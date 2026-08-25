r"""Tests for the post-Issue-#186 strict-1-arg shim.

The :class:`_BoundBoundaryOperator` shim pairs a realized 1-arg
:class:`LinearOperator` produced by
:class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer` with the LAW
it was realized from, and adds three surfaces:

* :meth:`apply(psi)` and :meth:`apply_transpose(psi)` — strict 1-arg
  passthroughs to the inner operator. The pre-#186 ``*_extra, **_kw``
  swallow affordance is gone; any test still calling
  ``bc.apply(psi, quad)`` would fail with ``TypeError``.
* the structural predicates and the function-space tags — delegate to
  the wrapped operator so the shim composes cleanly with other Wave-0
  primitives.
* :attr:`law`, and :attr:`kind` + ``__eq__`` against strings —
  :attr:`kind` preserves the legacy
  ``sn_mesh.bc["xmin"] == "reflective"`` comparison surface, now as a
  read-through of ``law``'s registry key.

**Fixture discipline (B2.0).** Every ``(inner, law)`` pair below is one
a realizer genuinely produces — vacuum with a mask / the zero map,
reflective with a permutation, albedo with a scaled identity. The shim
cannot check the pairing itself (realization is method-specific), so an
arbitrary pairing here would silently teach a false correspondence to
the next reader.

History
=======

The Wave-8/9 implementation carried an optional ``quadrature=``
kwarg that bound an ``AngularQuadrature`` and forwarded
``inner.apply(psi, bound_quad)`` to a legacy 2-arg
:class:`BoundaryTraceLaw` body. Issue #176 dropped that mode after
Issue #188 (curvilinear trace support) eliminated
the curvilinear-bypass code path. Issue #186 (B3 + β2) then dropped
the ``*_extra, **_kw`` swallow on :meth:`apply` since every remaining
caller is strict 1-arg. Campaign phase **B2.0** made ``law`` a required
constructor argument and turned ``kind`` from a stored string into a
property of that law.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry.boundary import (
    AlbedoBoundary,
    BoundaryTraceLaw,
    ReflectiveBoundary,
    VacuumInflow,
)
from orpheus.geometry.boundary._bound_compat import _BoundBoundaryOperator
from orpheus.numerics.operator import (
    IdentityOperator,
    PermutationOperator,
    ZeroOperator,
)
from tests._harness.references import mirror_partner_indices
from tests.sn._test_helpers import placeholder_materials


pytestmark = pytest.mark.l0


def test_apply_is_strict_1_arg():
    """``apply(psi)`` delegates to ``inner.apply(psi)``. Issue #186
    (B3 + β2) dropped the ``*_extra, **_kw`` swallow — extras now
    raise ``TypeError`` at the call site (matching the modernized
    1-arg :class:`LinearOperator` contract).
    """
    inner = IdentityOperator()
    shim = _BoundBoundaryOperator(inner, AlbedoBoundary(albedo=1.0))
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
    shim = _BoundBoundaryOperator(inner, ReflectiveBoundary())
    psi = np.arange(8.0).reshape(4, 2)
    np.testing.assert_array_equal(
        shim.apply_transpose(psi),
        inner.apply_transpose(psi),
    )
    # Strict 1-arg: extras raise.
    with pytest.raises(TypeError):
        shim.apply_transpose(psi, "extra", k=1)


def test_predicates_delegate_to_inner():
    """The structural predicates return the wrapped operator's truth —
    a permutation brings invertibility + adjointability; the zero map
    brings only adjointability (carve P4 rewire of the caps-delegation
    pin: same delegation contract, predicate spelling).
    """
    inner_perm = PermutationOperator(np.array([1, 0]), axis=0)
    perm_shim = _BoundBoundaryOperator(inner_perm, ReflectiveBoundary())
    assert perm_shim.is_invertible == inner_perm.is_invertible is True
    assert perm_shim.is_adjointable == inner_perm.is_adjointable is True

    zero_shim = _BoundBoundaryOperator(ZeroOperator(), VacuumInflow())
    assert zero_shim.is_invertible is False
    assert zero_shim.is_adjointable is True


def test_composes_with_operator_algebra():
    """The shim is a first-class :class:`LinearOperator` — it inherits
    the algebra dunders from :class:`LinearOperator`. ``shim + shim``
    builds an :class:`OperatorSum`; ``2.0 * shim`` builds a
    :class:`ScaledOperator`.
    """
    inner = IdentityOperator()
    shim = _BoundBoundaryOperator(inner, AlbedoBoundary(albedo=1.0))
    psi = np.ones((3, 2))
    # scalar multiplication via __mul__
    scaled = 2.0 * shim
    np.testing.assert_array_equal(scaled.apply(psi), 2.0 * psi)
    # sum via __add__ — note OperatorSum.apply forwards 1-arg only
    summed = shim + shim
    np.testing.assert_array_equal(summed.apply(psi), 2.0 * psi)


def test_kind_tag_supports_legacy_string_equality():
    """:attr:`kind` reads the LAW's registry key and ``__eq__`` compares
    it against strings — preserving the legacy SN-side
    ``sn_mesh.bc["xmin"] == "reflective"`` comparison that
    test_boundary_conditions.py + the BC-resolution diagnostic rely on.

    B2.0 turned ``kind`` from a constructor string into a property of
    :attr:`law`, so the tag can no longer be set to something the law
    does not say.
    """
    tagged = _BoundBoundaryOperator(ZeroOperator(), VacuumInflow())
    assert tagged == "vacuum"
    assert tagged != "reflective"

    # A law with no registry slot reports ``kind is None``, and the
    # string compare is then False against every string. Pre-B2.0 this
    # was reachable by omitting ``kind=``; it is now reachable only
    # through a law that never claimed a key — which is the honest
    # spelling of "this operator's law has no tag".
    class _Unregistered(BoundaryTraceLaw):  # no ``key=`` → not registered
        pass

    untagged = _BoundBoundaryOperator(IdentityOperator(), _Unregistered())
    assert untagged.kind is None
    assert (untagged == "vacuum") is False
    assert (untagged == "anything") is False


def test_kind_reads_the_registry_key_not_the_law_s_kind():
    r"""The read-through targets ``law.key``, NOT ``law.kind`` — measured.

    They agree for six of the seven laws and diverge for exactly one: a
    partially-reflecting :class:`ReflectiveBoundary` reports
    ``kind == "partial"`` (mirroring the ``BC("partial", albedo=…)``
    declaration vocabulary — the B0.1 ruling) while its ``key`` stays
    ``"reflective"`` for every albedo.

    Pre-B2.0 the shim stored ``law.key``, so **the key is the
    behaviour-preserving choice**; sourcing ``law.kind`` here would drop
    partially-reflecting faces out of
    ``sweep_schedule.reflective_faces``' ``== "reflective"`` set — a
    semantic change wearing a refactor's clothes. This leg is what
    reddens if someone "tidies" the property to the more obvious name.
    """
    partial = ReflectiveBoundary(albedo=0.7)
    assert partial.kind == "partial"          # the LAW's own answer
    assert type(partial).key == "reflective"  # the REGISTRY's answer

    shim = _BoundBoundaryOperator(
        PermutationOperator(np.array([1, 0]), axis=0), partial,
    )
    assert shim.kind == "reflective"
    assert shim == "reflective"
    assert shim != "partial"


def test_shim_carries_the_law_it_was_realized_from():
    """B2.0's whole content: the descriptor survives realization.

    Before B2.0 ``sn_mesh.bc[face]`` was a realized operator plus a
    string, so a consumer could ask what the face was *declared* as but
    not what its law *does*. The five production string-dispatch sites
    are that gap. The law's two affine factors are reachable here.
    """
    law = ReflectiveBoundary(axis="y", albedo=0.7)
    shim = _BoundBoundaryOperator(
        PermutationOperator(np.array([1, 0]), axis=0), law,
    )
    assert shim.law is law
    # The structural questions the string sites are really asking.
    assert shim.law.geometry_map.permutes_ordinates is True
    assert shim.law.response_kernel.amplitude == 0.7


def test_space_tags_forward_to_inner():
    """``domain`` / ``codomain`` delegate like every other wrapper.

    A transparent handle that reported the base ``None`` default would
    make composition checks *skip* silently rather than run (see
    ``LinearOperator.domain``). The shim forwarded eight members and not
    these two until B2.0.
    """
    inner = PermutationOperator(np.array([1, 0]), axis=0)
    shim = _BoundBoundaryOperator(inner, ReflectiveBoundary())
    assert shim.domain is getattr(inner, "domain", None)
    assert shim.codomain is getattr(inner, "codomain", None)


def test_shim_remains_hashable():
    """``__hash__`` falls back to identity so the shim is usable as a
    dict key / set element. Two shims wrapping the same inner are
    distinct (different id), preserving the standard Python identity-
    of-instance semantics.
    """
    inner = ZeroOperator()
    a = _BoundBoundaryOperator(inner, VacuumInflow())
    b = _BoundBoundaryOperator(inner, VacuumInflow())
    assert hash(a) != hash(b)  # distinct instances
    d = {a: 1, b: 2}
    assert d[a] == 1
    assert d[b] == 2


def test_shim_is_not_re_exported_from_package():
    """The shim is internal to the SN-side wiring — it MUST NOT appear
    in :mod:`orpheus.geometry.boundary`'s public surface, since no
    consumer outside ``SNMesh.realize_boundary_law`` has a legitimate
    reason to wrap operators in it. Pinning the lack of re-export
    prevents accidental promotion to public API.
    """
    import orpheus.geometry.boundary as boundary_pkg

    assert "_BoundBoundaryOperator" not in getattr(boundary_pkg, "__all__", [])
    assert not hasattr(boundary_pkg, "_BoundBoundaryOperator")


def test_constructor_signature_is_inner_plus_law():
    """The signature pin: ``__init__(inner, law)`` — nothing else.

    Two retirements are pinned at once.

    * Issue #176 / C176.1: the Wave-8/9 dual-mode ``quadrature=`` kwarg
      and its ``_quadrature`` attribute are gone.
    * Campaign B2.0: ``kind=`` is gone. A realized boundary law that
      cannot say which law it realizes is the state B2 exists to
      delete, so ``law`` is required and the tag is derived — a
      caller can no longer hand the shim a string that contradicts it.
    """
    inner = ZeroOperator()
    shim = _BoundBoundaryOperator(inner, VacuumInflow())
    assert not hasattr(shim, "_quadrature")
    with pytest.raises(TypeError):
        _BoundBoundaryOperator(inner, quadrature="not_accepted")  # type: ignore[call-arg]
    with pytest.raises(TypeError):
        _BoundBoundaryOperator(inner, kind="vacuum")  # type: ignore[call-arg]
    with pytest.raises(TypeError):
        _BoundBoundaryOperator(inner)  # type: ignore[call-arg]


# ═══════════════════════════════════════════════════════════════════════
# Issue #188 / C188.3 — curvilinear realizer wiring
# ═══════════════════════════════════════════════════════════════════════
#
# These tests pin the SNMesh-side production behaviour: a 1-D curvilinear
# mesh routes its BCs through SNBoundaryRealizer just like Cartesian, and
# the 1-D y-face placeholders are realized through SNMethodSpace.minimal.


class Test188WiringContracts:
    """Curvilinear ``SNMesh.realize_boundary_law`` + 1-D y-placeholder
    contracts (the per-face arm was named ``_resolve_one`` until the
    #290 P7b ``TransportMethod`` carve)."""

    def test_curvilinear_realize_boundary_law_routes_through_realizer(self):
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
        from orpheus.numerics.operator import ZeroMorphism
        from orpheus.sn.mesh.augmented_mesh import SNMesh
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
        # primitive.  RE-POSED at campaign phase B3.2: with the law's domain
        # narrowed to Γ₊, vacuum's whole content (R = 0) IS the zero map
        # Γ₊ → Γ₋, so there is no full-face projector left to decompose into
        # the §16A.10 ``B = G_patch ⊗ K_omega ⊗ K_g`` tensor-product form.
        # The pre-B3.2 pin drilled into that TP for an
        # ``IncomingOrdinateMaskTensor``; the claim this test exists for —
        # curvilinear meshes route through the SAME realizer — is unchanged,
        # only the object at the end of the route moved.
        assert isinstance(sn.bc["xmax"], _BoundBoundaryOperator)
        assert isinstance(sn.bc["xmax"].inner, ZeroMorphism)
        # What the resolution MEANS — the law that landed on the face.
        assert isinstance(sn.bc["xmax"].law, VacuumInflow)
        # …and the legacy string-comparison surface still answers for it,
        # on a REAL resolved mesh (the hand-built-shim leg above cannot
        # show that the production path still feeds it).
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

        ⭐ RE-POSED at **B3.2**, and this fixture is one of the TWO
        discriminating fixtures for the narrowing's new index remap. The
        realized permutation is now on the reduced ordinate axis, i.e.
        ``local_positions(mirror_partner[inflow], outflow)``, and on a
        SLAB the mirror REVERSES order: with ``gauss_legendre(4)`` at ``xmax``,
        ``perm[inflow] = [3, 2]`` maps to local ``[1, 0]`` while the naive
        ``arange(2)`` would give ``[0, 1]``. So the slab — not a curvilinear or
        2-D mesh — is where a hand-written ``arange`` in the REFLECTIVE
        narrowing bites. (Its sibling trap, the SCHEDULE-SPLIT remap, is the
        opposite: correct in 1-D, wrong in 2-D — gated by RG-5 in
        ``tests/sn/operators/test_sn_boundary_operator.py``. Neither fixture
        covers the other; both gates are load-bearing.)
        """
        from orpheus.geometry import BC, CoordSystem, Mesh1D
        from orpheus.numerics.operator import TensorProductOperator
        from orpheus.sn.mesh.augmented_mesh import SNMesh
        from orpheus.numerics.quadrature import Quadrature
        from tests.sn._test_helpers import local_positions

        mesh = Mesh1D(
            edges=np.linspace(0.0, 1.0, 5),
            mat_ids=np.zeros(4, dtype=int),
            coord=CoordSystem.CARTESIAN,
            bc_left=BC("reflective"),
            bc_right=BC("reflective"),
        )
        quad = Quadrature.gauss_legendre(4)
        sn = SNMesh(mesh, quad, placeholder_materials())

        discriminating = 0
        for face in ("xmin", "xmax"):
            assert isinstance(sn.bc[face], _BoundBoundaryOperator), face
            assert isinstance(sn.bc[face].inner, TensorProductOperator), face
            angular = sn.bc[face].inner.ops[0]
            assert isinstance(angular, PermutationOperator), face
            inflow = sn.angular_trace.inflow_indices_for_face(face)
            outflow = sn.angular_trace.outflow_indices_for_face(face)
            # GL1D x-reflection pairs mu with -mu — the flip permutation,
            # narrowed to the face's half-traces.
            expected = local_positions(
                mirror_partner_indices(quad, "x")[inflow], outflow,
            )
            np.testing.assert_array_equal(angular.perm, expected)
            if not np.array_equal(expected, np.arange(expected.size)):
                discriminating += 1
            # The law that landed, and the legacy tag surface it feeds.
            assert isinstance(sn.bc[face].law, ReflectiveBoundary)
            assert sn.bc[face] == "reflective"
        # ACTIVATION guard: without it this gate could silently decay into a
        # fixture where the naive ``arange`` happens to be right, and stop
        # testing the remap at all (vv Mode 8, class 7 — a decayed gate).
        if discriminating == 0:
            pytest.fail(
                "no slab face has a reduced permutation that differs from "
                "arange — this fixture no longer discriminates the reflective "
                "narrowing's to_local remap."
            )
