r"""Tests for the Wave-8 transitional 2-arg shim (C8.2).

The :class:`_BoundBoundaryOperator` shim wraps a realized 1-arg
:class:`LinearOperator` so the 13 production call sites at
:mod:`orpheus.sn.sweep` + :mod:`orpheus.sn.operator` (which still call
``bc.apply(psi, quad)`` with 2 args) keep working through Wave 9.

These tests pin:

* ``apply(psi, *anything, **anything)`` forwards to ``inner.apply(psi)``
  — extra args are swallowed.
* ``apply_transpose`` symmetrical for inner ops that support it.
* :attr:`capabilities` property delegates to the wrapped operator.
* The shim composes cleanly with the operator-algebra dunders inherited
  from :class:`LinearOperatorMixin` — proves it is a first-class
  :class:`LinearOperator`.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry.boundary._bound_compat import _BoundBoundaryOperator
from orpheus.numerics.operator import (
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    IdentityOperator,
    PermutationOperator,
    ZeroOperator,
)


pytestmark = pytest.mark.l0


def test_apply_forwards_and_swallows_extra_args():
    """``apply(psi, *_extra, **_kw)`` delegates to ``inner.apply(psi)`` and
    ignores any positional / keyword extras (typically a quadrature
    argument passed by the legacy 2-arg call sites).
    """
    inner = IdentityOperator()
    shim = _BoundBoundaryOperator(inner)
    psi = np.arange(12.0).reshape(4, 3)
    # 1-arg
    np.testing.assert_array_equal(shim.apply(psi), psi)
    # 2-arg — extra "quadrature-like" sentinel is swallowed
    sentinel = object()
    np.testing.assert_array_equal(shim.apply(psi, sentinel), psi)
    # kwargs swallowed too
    np.testing.assert_array_equal(
        shim.apply(psi, quad=sentinel, extra=42), psi,
    )


def test_apply_transpose_forwards_when_inner_supports_it():
    """For an inner operator with ``apply_transpose`` (e.g. a
    permutation), the shim's :meth:`apply_transpose` returns the
    same result as ``inner.apply_transpose``.
    """
    perm = np.array([2, 0, 1, 3])
    inner = PermutationOperator(perm, axis=0)
    shim = _BoundBoundaryOperator(inner)
    psi = np.arange(8.0).reshape(4, 2)
    np.testing.assert_array_equal(
        shim.apply_transpose(psi),
        inner.apply_transpose(psi),
    )
    # And again ignoring extra args
    np.testing.assert_array_equal(
        shim.apply_transpose(psi, "ignored", k=1),
        inner.apply_transpose(psi),
    )


def test_capabilities_delegate_to_inner():
    """:attr:`capabilities` returns whatever the wrapped operator
    advertises — a permutation brings both apply + apply_transpose; a
    mask-only operator brings only apply.
    """
    perm_shim = _BoundBoundaryOperator(PermutationOperator(np.array([1, 0]), axis=0))
    assert perm_shim.capabilities == frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})

    zero_shim = _BoundBoundaryOperator(ZeroOperator())
    # ZeroOperator advertises {CAP_APPLY, CAP_APPLY_TRANSPOSE}
    assert CAP_APPLY in zero_shim.capabilities


def test_composes_with_operator_algebra():
    """The shim is a first-class :class:`LinearOperator` — it inherits
    the algebra dunders from :class:`LinearOperatorMixin`. ``shim + shim``
    builds an :class:`OperatorSum`; ``2.0 * shim`` builds a
    :class:`ScaledOperator`. This is the load-bearing piece for the
    realizer's mixed-BC recursive path (when MixedBoundaryOperator is
    composed with a wrapped realized op).
    """
    inner = IdentityOperator()
    shim = _BoundBoundaryOperator(inner)
    psi = np.ones((3, 2))
    # scalar multiplication via __mul__
    scaled = 2.0 * shim
    np.testing.assert_array_equal(scaled.apply(psi), 2.0 * psi)
    # sum via __add__ — note OperatorSum.apply forwards 1-arg only
    summed = shim + shim
    np.testing.assert_array_equal(summed.apply(psi), 2.0 * psi)


def test_kind_tag_supports_legacy_string_equality():
    """The shim accepts an optional ``kind`` tag and implements
    ``__eq__`` against strings — preserves the legacy SN-side
    ``sn_mesh.bc_xmin == "reflective"`` comparison that test_boundary_conditions.py
    + the BC-resolution diagnostic rely on. Without a kind tag the
    comparison returns ``NotImplemented`` (so ``shim == "x"`` is False).
    """
    inner = IdentityOperator()
    tagged = _BoundBoundaryOperator(inner, kind="vacuum")
    assert tagged == "vacuum"
    assert tagged != "reflective"

    # Untagged shim: string compare is NotImplemented → False from ==
    untagged = _BoundBoundaryOperator(inner)
    assert (untagged == "vacuum") is False
    assert (untagged == "anything") is False


def test_shim_remains_hashable():
    """``__hash__`` falls back to identity so the shim is usable as a
    dict key / set element. Two shims wrapping the same inner are
    distinct (different id), preserving the standard Python identity-
    of-instance semantics.
    """
    inner = IdentityOperator()
    a = _BoundBoundaryOperator(inner, kind="vacuum")
    b = _BoundBoundaryOperator(inner, kind="vacuum")
    assert hash(a) != hash(b)  # distinct instances
    d = {a: 1, b: 2}
    assert d[a] == 1
    assert d[b] == 2


def test_shim_is_not_re_exported_from_package():
    """The shim is internal to the SN-side wiring — it MUST NOT appear
    in :mod:`orpheus.geometry.boundary`'s public surface, since no
    consumer outside ``SNMesh._resolve_bcs`` has a legitimate reason
    to wrap operators in it. Pinning the lack of re-export prevents
    accidental promotion to public API.
    """
    import orpheus.geometry.boundary as boundary_pkg

    assert "_BoundBoundaryOperator" not in getattr(boundary_pkg, "__all__", [])
    assert not hasattr(boundary_pkg, "_BoundBoundaryOperator")


# ═══════════════════════════════════════════════════════════════════════
# Wave 9 (C9.0) — bound-quadrature mode for the curvilinear path
# ═══════════════════════════════════════════════════════════════════════
#
# When the shim wraps a LEGACY 2-arg :class:`BoundaryTraceLaw` instance
# (curvilinear path in :meth:`SNMesh._resolve_one` — Wave 2 has not yet
# implemented :class:`InflowTraceSpace` for spherical / cylindrical
# meshes), the shim is built with ``quadrature=<bound quad>`` so that
# the uniform 1-arg ``bc.apply(psi)`` contract holds across ALL ``bc_*``
# attributes.  The shim then forwards
# ``inner.apply(psi, self._quadrature)`` — bit-identical to the
# pre-Wave-9 ``bc.apply(psi, quad)`` direct call.


class _Recording2ArgInner:
    """Inner stub that records the (psi, quadrature) it was called with.

    Mimics the legacy ``BoundaryTraceLaw.apply(psi, quadrature)`` 2-arg
    signature so we can pin that the shim's bound-quadrature forwarding
    delivers exactly what the legacy call sites used to pass.
    """

    def __init__(self) -> None:
        self.apply_calls: list[tuple[object, object]] = []
        self.transpose_calls: list[tuple[object, object]] = []
        self.capabilities = frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})

    def apply(self, psi, quadrature):
        self.apply_calls.append((psi, quadrature))
        return psi  # echo for assertions

    def apply_transpose(self, psi, quadrature):
        self.transpose_calls.append((psi, quadrature))
        return psi


class Test2ArgForwarding:
    """Bound-quadrature mode for the curvilinear ``SNMesh._resolve_one`` path."""

    def test_apply_forwards_bound_quadrature_to_inner(self):
        """``apply(psi)`` forwards ``inner.apply(psi, self._quadrature)``
        when the shim was built with a bound quadrature. Pins that the
        bit-identical legacy call ``bc.apply(psi, quad)`` survives the
        Wave 9 1-arg migration.
        """
        inner = _Recording2ArgInner()
        sentinel_quad = object()  # identity-checked below
        shim = _BoundBoundaryOperator(
            inner, quadrature=sentinel_quad, kind="reflective",
        )
        psi = np.arange(6.0).reshape(3, 2)

        # 1-arg call (Wave 9 target signature)
        np.testing.assert_array_equal(shim.apply(psi), psi)
        assert len(inner.apply_calls) == 1
        psi_recv, quad_recv = inner.apply_calls[0]
        np.testing.assert_array_equal(psi_recv, psi)
        # Identity, not just equality — pins that the bound quadrature
        # is the SAME object the realizer passed at construction time.
        assert quad_recv is sentinel_quad

    def test_apply_swallows_extra_legacy_quad_when_bound(self):
        """During the mid-Wave-9 transition some call sites may still
        pass ``quad`` as a 2nd positional; the shim swallows it and
        STILL forwards its own bound quadrature, NOT the one the caller
        passed. This keeps the bound-quadrature contract authoritative.
        """
        inner = _Recording2ArgInner()
        bound_quad = object()
        stale_quad = object()  # what a not-yet-migrated caller passes
        shim = _BoundBoundaryOperator(inner, quadrature=bound_quad)
        psi = np.ones((2, 2))

        shim.apply(psi, stale_quad)
        _, quad_recv = inner.apply_calls[-1]
        assert quad_recv is bound_quad
        assert quad_recv is not stale_quad

        # kwargs swallowed too
        shim.apply(psi, quad=stale_quad, extra=42)
        _, quad_recv = inner.apply_calls[-1]
        assert quad_recv is bound_quad

    def test_apply_transpose_forwards_bound_quadrature(self):
        """Symmetric to ``apply``: ``apply_transpose(psi)`` forwards
        ``inner.apply_transpose(psi, self._quadrature)``.
        """
        inner = _Recording2ArgInner()
        bound_quad = object()
        shim = _BoundBoundaryOperator(inner, quadrature=bound_quad)
        psi = np.arange(4.0).reshape(2, 2)

        np.testing.assert_array_equal(shim.apply_transpose(psi), psi)
        assert len(inner.transpose_calls) == 1
        psi_recv, quad_recv = inner.transpose_calls[0]
        np.testing.assert_array_equal(psi_recv, psi)
        assert quad_recv is bound_quad

    def test_quadrature_none_preserves_wave8_1arg_path(self):
        """Default ``quadrature=None`` still forwards 1-arg
        ``inner.apply(psi)`` — the Wave 8 Cartesian-path behaviour
        MUST remain bit-identical. This is the regression contract
        against the 7 pre-existing tests above.
        """
        inner = IdentityOperator()
        shim = _BoundBoundaryOperator(inner)  # no quadrature kwarg
        assert shim._quadrature is None
        psi = np.arange(6.0).reshape(2, 3)
        # Bit-identical to inner.apply(psi)
        np.testing.assert_array_equal(shim.apply(psi), inner.apply(psi))

    def test_curvilinear_resolve_one_routes_through_realizer(self):
        """Issue #188 / C188.3: a 1-D spherical :class:`SNMesh`
        exposes ``bc_left`` / ``bc_right`` as
        :class:`_BoundBoundaryOperator` shims wrapping REALIZED
        1-arg operators (no bound quadrature). C188.1+C188.2
        extended :class:`InflowTraceSpace` to all 1-D coord
        systems, so the realizer-then-shim path is now used
        uniformly across Cartesian / spherical / cylindrical
        meshes.

        Compatibility against the pre-C188.3 bound-quadrature path
        is verified end-to-end by the curvilinear SN regression
        suite (``tests/sn/test_spherical.py`` etc.); here we pin
        the structural contract.
        """
        from orpheus.geometry import BC, CoordSystem, Mesh1D
        from orpheus.numerics.operator import (
            IncomingOrdinateMaskTensor,
            PermutationOperator as _PO,
        )
        from orpheus.sn.geometry import SNMesh
        from orpheus.sn.quadrature import GaussLegendre1D

        mesh = Mesh1D(
            edges=np.linspace(0.0, 1.0, 5),
            mat_ids=np.zeros(4, dtype=int),
            coord=CoordSystem.SPHERICAL,
            bc_left=BC("reflective"),
            bc_right=BC("vacuum"),
        )
        quad = GaussLegendre1D.create(8)
        sn = SNMesh(mesh, quad)

        # Curvilinear now builds traces (C188.1+C188.2 lifted guard).
        assert sn._inflow_trace is not None
        # Realizer-path shims: no bound quadrature.
        assert isinstance(sn.bc_left, _BoundBoundaryOperator)
        assert isinstance(sn.bc_right, _BoundBoundaryOperator)
        assert sn.bc_left._quadrature is None
        assert sn.bc_right._quadrature is None
        # Inner is the realized primitive: PermutationOperator for
        # reflective; IncomingOrdinateMaskTensor for vacuum.
        assert isinstance(sn.bc_left.inner, _PO)
        assert isinstance(sn.bc_right.inner, IncomingOrdinateMaskTensor)
        # Kind tags survive for the legacy string-comparison surface.
        assert sn.bc_left == "reflective"
        assert sn.bc_right == "vacuum"

    def test_1d_y_face_placeholders_realized_through_minimal_space(self):
        """Issue #188 / C188.3: the 1-D ``bc_ymin`` / ``bc_ymax``
        placeholders route through :class:`SNBoundaryRealizer` with
        :meth:`SNMethodSpace.minimal` — no bound quadrature is
        needed, since the realized op is already 1-arg. The 1-D
        trace space cannot service
        ``inflow_indices_for_face('ymin')`` (its face_names are
        ``("left", "right")`` only) but the realizer's
        :class:`ReflectiveBoundary` branch does NOT consume
        inflow_indices, only ``law.axis`` and
        ``quad.reflection_index``. For :class:`GaussLegendre1D`,
        ``reflection_index("y")`` returns the identity permutation
        (every ordinate is its own partner because ``mu_y == 0``),
        so the realized op is a no-op
        :class:`PermutationOperator`.
        """
        from orpheus.geometry import BC, CoordSystem, Mesh1D
        from orpheus.sn.geometry import SNMesh
        from orpheus.sn.quadrature import GaussLegendre1D

        mesh = Mesh1D(
            edges=np.linspace(0.0, 1.0, 5),
            mat_ids=np.zeros(4, dtype=int),
            coord=CoordSystem.CARTESIAN,
            bc_left=BC("vacuum"),
            bc_right=BC("vacuum"),
        )
        quad = GaussLegendre1D.create(4)
        sn = SNMesh(mesh, quad)

        # y-face placeholders shim-wrap a REALIZED PermutationOperator.
        assert isinstance(sn.bc_ymin, _BoundBoundaryOperator)
        assert isinstance(sn.bc_ymax, _BoundBoundaryOperator)
        assert isinstance(sn.bc_ymin.inner, PermutationOperator)
        assert isinstance(sn.bc_ymax.inner, PermutationOperator)
        # The realized op is 1-arg; no bound quadrature on the shim.
        assert sn.bc_ymin._quadrature is None
        assert sn.bc_ymax._quadrature is None
        # Kind preserved for the legacy string-compare surface.
        assert sn.bc_ymin == "reflective"
        assert sn.bc_ymax == "reflective"
