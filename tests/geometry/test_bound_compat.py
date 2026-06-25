r"""Tests for the post-Issue-#186 strict-1-arg shim.

The :class:`_BoundBoundaryOperator` shim wraps a realized 1-arg
:class:`LinearOperator` produced by
:class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer` and adds
two thin surfaces:

* :meth:`apply(psi)` and :meth:`apply_transpose(psi)` — strict 1-arg
  passthroughs to the inner operator. The pre-#186 ``*_extra, **_kw``
  swallow affordance is gone; any test still calling
  ``bc.apply(psi, quad)`` would fail with ``TypeError``.
* :attr:`capabilities` — delegates to the wrapped operator so the
  shim composes cleanly with other Wave-0 primitives.
* :attr:`kind` + ``__eq__`` against strings — preserves the legacy
  ``sn_mesh.bc["xmin"] == "reflective"`` comparison surface.

History
=======

The Wave-8/9 implementation carried an optional ``quadrature=``
kwarg that bound an :class:`AngularQuadrature` and forwarded
``inner.apply(psi, bound_quad)`` to a legacy 2-arg
:class:`BoundaryTraceLaw` body. Issue #176 dropped that mode after
Issue #188 (curvilinear trace support) eliminated
the curvilinear-bypass code path. Issue #186 (B3 + β2) then dropped
the ``*_extra, **_kw`` swallow on :meth:`apply` since every remaining
caller is strict 1-arg.
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
from tests.sn._test_helpers import placeholder_materials


pytestmark = pytest.mark.l0


def test_apply_is_strict_1_arg():
    """``apply(psi)`` delegates to ``inner.apply(psi)``. Issue #186
    (B3 + β2) dropped the ``*_extra, **_kw`` swallow — extras now
    raise ``TypeError`` at the call site (matching the modernized
    1-arg :class:`LinearOperator` contract).
    """
    inner = IdentityOperator()
    shim = _BoundBoundaryOperator(inner)
    psi = np.arange(12.0).reshape(4, 3)
    # 1-arg works
    np.testing.assert_array_equal(shim.apply(psi), psi)
    # 2-arg raises
    with pytest.raises(TypeError):
        shim.apply(psi, object())
    # extra kwargs raise
    with pytest.raises(TypeError):
        shim.apply(psi, quad=object())


def test_apply_transpose_forwards_when_inner_supports_it():
    """For an inner operator with ``apply_transpose`` (e.g. a
    permutation), the shim's :meth:`apply_transpose` returns the
    same result as ``inner.apply_transpose``. Issue #186 cleanup:
    strict 1-arg; extras raise ``TypeError``.
    """
    perm = np.array([2, 0, 1, 3])
    inner = PermutationOperator(perm, axis=0)
    shim = _BoundBoundaryOperator(inner)
    psi = np.arange(8.0).reshape(4, 2)
    np.testing.assert_array_equal(
        shim.apply_transpose(psi),
        inner.apply_transpose(psi),
    )
    # Strict 1-arg: extras raise.
    with pytest.raises(TypeError):
        shim.apply_transpose(psi, "extra", k=1)


def test_capabilities_delegate_to_inner():
    """:attr:`capabilities` returns whatever the wrapped operator
    advertises — a permutation brings apply + apply_transpose + solve
    (post Issue #150 closeout); a mask-only operator brings only apply.

    The shim's contract is delegation; the exact inner capability set
    is the inner operator's concern.
    """
    inner_perm = PermutationOperator(np.array([1, 0]), axis=0)
    perm_shim = _BoundBoundaryOperator(inner_perm)
    # Delegation: shim advertises exactly what inner advertises.
    assert perm_shim.capabilities == inner_perm.capabilities

    zero_shim = _BoundBoundaryOperator(ZeroOperator())
    # ZeroOperator advertises {CAP_APPLY, CAP_APPLY_TRANSPOSE}
    assert CAP_APPLY in zero_shim.capabilities


def test_composes_with_operator_algebra():
    """The shim is a first-class :class:`LinearOperator` — it inherits
    the algebra dunders from :class:`LinearOperator`. ``shim + shim``
    builds an :class:`OperatorSum`; ``2.0 * shim`` builds a
    :class:`ScaledOperator`.
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
    ``sn_mesh.bc["xmin"] == "reflective"`` comparison that
    test_boundary_conditions.py + the BC-resolution diagnostic rely
    on. Without a kind tag the comparison returns ``NotImplemented``
    (so ``shim == "x"`` is False).
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


def test_shim_has_no_quadrature_attribute_after_176():
    """Issue #176 / C176.1: the Wave-8/9 dual-mode shim's
    ``quadrature=`` kwarg and ``_quadrature`` attribute are gone.
    Pin the post-cleanup signature: ``__init__(inner, kind=None)``
    has no quadrature parameter, and instances carry no
    ``_quadrature`` attribute.
    """
    inner = IdentityOperator()
    shim = _BoundBoundaryOperator(inner, kind="vacuum")
    assert not hasattr(shim, "_quadrature")
    with pytest.raises(TypeError):
        _BoundBoundaryOperator(inner, quadrature="not_accepted")  # type: ignore[call-arg]


# ═══════════════════════════════════════════════════════════════════════
# Issue #188 / C188.3 — curvilinear realizer wiring
# ═══════════════════════════════════════════════════════════════════════
#
# These tests pin the SNMesh-side production behaviour: a 1-D curvilinear
# mesh routes its BCs through SNBoundaryRealizer just like Cartesian, and
# the 1-D y-face placeholders are realized through SNMethodSpace.minimal.


class Test188WiringContracts:
    """Curvilinear ``SNMesh._resolve_one`` + 1-D y-placeholder contracts."""

    def test_curvilinear_resolve_one_routes_through_realizer(self):
        """A 1-D spherical :class:`SNMesh` has exactly ONE boundary —
        the outer radius (``xmax`` / ``bc_right``) — realized through
        :class:`SNBoundaryRealizer` to a 1-arg operator. The pole r=0
        is the angular closure's regularity condition, NOT a BC face,
        so ``bc_left`` / ``bc_xmin`` are ``None`` (the mesh's
        ``bc_left`` declaration is moot — the centreline is always
        symmetric by geometry).

        Compatibility against the pre-C188.3 bound-quadrature path
        is verified end-to-end by the curvilinear SN regression
        suite (``tests/sn/test_spherical.py`` etc.); here we pin
        the structural contract.
        """
        from orpheus.geometry import BC, CoordSystem, Mesh1D
        from orpheus.numerics.operator import (
            IncomingOrdinateMaskTensor,
            TensorProductOperator,
        )
        from orpheus.sn.geometry import SNMesh
        from orpheus.numerics.quadrature import Quadrature

        mesh = Mesh1D(
            edges=np.linspace(0.0, 1.0, 5),
            mat_ids=np.zeros(4, dtype=int),
            coord=CoordSystem.SPHERICAL,
            bc_left=BC("reflective"),
            bc_right=BC("vacuum"),
        )
        quad = Quadrature.gauss_legendre(8)
        sn = SNMesh(mesh, quad, placeholder_materials())

        # ONE boundary: the unified trace carries only the outer face.
        assert sn._trace is not None
        assert sn._trace.face_names == ("xmax",)
        # No inner-face entry at the pole (C4 / #220 — structurally
        # absent from the face-name-keyed ``bc`` dict, not ``None``).
        assert set(sn.bc) == {"xmax"}
        # Outer face: realizer-path shim wrapping the realized vacuum
        # primitive.  Post-Wave-T.1 the realizer lifts the vacuum mask
        # into the streaming tensor-product form
        # ``IncomingOrdinateMaskTensor(axis=0) ⊗ IdentityOperator`` (the
        # §16A.10 ``B = G_patch ⊗ K_omega ⊗ K_g`` decomposition, where
        # the only non-trivial factor is the ordinate-axis mask).  Drill
        # into the tensor product to pin that the realized angular factor
        # is STILL the vacuum mask — the structure changed, the semantics
        # did not.
        assert isinstance(sn.bc["xmax"], _BoundBoundaryOperator)
        assert isinstance(sn.bc["xmax"].inner, TensorProductOperator)
        ordinate_factor = sn.bc["xmax"].inner.ops[0]
        assert isinstance(ordinate_factor, IncomingOrdinateMaskTensor)
        # Kind tag survives for the legacy string-comparison surface.
        assert sn.bc["xmax"] == "vacuum"

    def test_1d_reflective_faces_realized_as_permutation_tp(self):
        """1-D slab reflective faces shim-wrap a REALIZED reflective
        primitive. Post-D-B+1 the realizer lifts the specular
        permutation into the streaming tensor-product form
        ``PermutationOperator(axis=0) ⊗ IdentityOperator`` (the
        §16A.10 ``B = G_patch ⊗ K_omega ⊗ K_g`` decomposition;
        albedo == 1.0 so the bare-permutation TP fast path is
        taken). Drill into the tensor product to pin that the
        realized angular factor is STILL the reflection permutation
        — the structure changed, the semantics did not.

        HISTORY: this pin originally lived on the degenerate 1-D
        ``bc_ymin`` / ``bc_ymax`` placeholders (Issue #188 /
        C188.3), which C4 (#220) retired — a 1-D mesh's ``bc``
        inventory now carries NO y-entries (no production code ever
        read the placeholders). The structural claim — realized
        reflective law in TP form behind the shim — is geometry
        content, so it moved onto the slab's REAL ``xmin`` /
        ``xmax`` faces.
        """
        from orpheus.geometry import BC, CoordSystem, Mesh1D
        from orpheus.numerics.operator import TensorProductOperator
        from orpheus.sn.geometry import SNMesh
        from orpheus.numerics.quadrature import Quadrature

        mesh = Mesh1D(
            edges=np.linspace(0.0, 1.0, 5),
            mat_ids=np.zeros(4, dtype=int),
            coord=CoordSystem.CARTESIAN,
            bc_left=BC("reflective"),
            bc_right=BC("reflective"),
        )
        quad = Quadrature.gauss_legendre(4)
        sn = SNMesh(mesh, quad, placeholder_materials())

        for face in ("xmin", "xmax"):
            assert isinstance(sn.bc[face], _BoundBoundaryOperator), face
            assert isinstance(sn.bc[face].inner, TensorProductOperator), face
            angular = sn.bc[face].inner.ops[0]
            assert isinstance(angular, PermutationOperator), face
            # GL1D x-reflection pairs mu with -mu — the flip permutation.
            np.testing.assert_array_equal(
                angular.perm, quad.reflection_index("x"),
            )
            # Kind preserved for the legacy string-compare surface.
            assert sn.bc[face] == "reflective"
