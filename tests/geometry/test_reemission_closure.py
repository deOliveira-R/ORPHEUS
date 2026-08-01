r"""Campaign phase **B3.4b** — the re-emission closure, and the discipline it enforces.

The user's 2026-08-01 ruling (`.claude/plans/b3_domain_narrowing_crosswalk.md`
§11.1) put the specular pairing of an absorbing *surface* in :math:`R`, not
:math:`G`: ``AlbedoBoundary.geometry_map`` is ``IdentityMap()``
**unconditionally**, and a :class:`ReemissionClosure` selects the response
kernel. This module carries the four claims that ruling creates and that have no
older home:

1. **The ≡ theorems** (:class:`TestEquivalenceTheorems`) —
   ``albedo(α, SpecularReturn(a)) ≡ reflective(a, α)`` and
   ``albedo(α, IsotropicReturn(a, s)) ≡ white(a, s, α)`` as realized operators.
2. **The independent-expression anchors**
   (:class:`TestSpecularAgainstAnIndependentExpression`,
   :class:`TestDiffuseAgainstAnIndependentExpression`) — the gates that actually
   catch a bug, because (1) cannot.
3. **The exactly-one-of-G-R invariant** (:class:`TestExactlyOneFactorIsNonTrivial`).
4. **The refusal, both directions** (:class:`TestBareAlbedoIsRefusedBySN`,
   :class:`TestDiffusionIsUnmovedAndClosureBlind`), plus the α=0 zero map the
   carve opened (:class:`TestZeroAmplitudeRealizesTheNarrowedZeroMap`).

Why (1) is NOT the verification, and (2) is
===========================================

The §12.2 design argues the ≡ theorems "hold because the two routes literally
execute the same body". That is true, and it is exactly why the equivalence gate
is **structurally non-independent**: a shared body carries a shared bug, so
:class:`TestEquivalenceTheorems` is GREEN under any error *inside* the body. It
catches ARGUMENT drift (a wrong axis, a dropped α, a wrong ``law_key``) and
nothing else — necessary, never sufficient (``vv-principles`` §1).

The catchers are the two anchor classes, whose references are written from
``quadrature.reflection_index`` / ``quadrature.weights`` × ``omega_dot_n`` and
import nothing from :mod:`orpheus.sn.boundary.realizer` above the ``np.take`` /
``np.tensordot`` line.

⭐ The two routes need COMPLEMENTARY fixtures — neither list covers the other
============================================================================

MEASURED 2026-08-01. The specular route's canonical mutation (the ``to_local``
remap written as a bare ``arange``) and the diffuse route's canonical mutation
(the B3.4a ``> 0.0`` classification twin) are discriminated by **disjoint**
quadratures:

.. list-table::
   :header-rows: 1

   * - mutation class
     - DISCRIMINATES on
     - BLIND on
   * - specular remap (``to_local`` → ``arange``)
     - ``gauss_legendre(4)``, ``gauss_legendre(8)``, ``lebedev(17)``
     - ``product(2,4)``, ``level_symmetric(6)``
   * - diffuse classification (``TANGENTIAL_EPS`` → ``> 0.0``)
     - ``product(2,4)``
     - ``gauss_legendre``, ``lebedev(17)``, ``level_symmetric(6)``

On ``product(2,4)`` and ``level_symmetric(6)`` the *positional* pairing a bare
albedo used to perform happens to EQUAL the specular pairing
(``perm[sorted(inflow)] == sorted(outflow)``); on the slab quadratures the mirror
REVERSES order and on ``lebedev(17)`` it reorders four entries. So both routes are
parametrised over the FULL fixture list and the blind rows are kept as documented
controls — pruning either one is designed-green. This is the crosswalk §9.1
finding ("``to_local``, two sites, complementary fixtures") reproduced from the
albedo side.

V&V tags
--------

``@pytest.mark.foundation`` throughout: these are SOFTWARE invariants about
operator construction and factor typing, with no theory-page ``:label:`` to
point at. No ``verifies(...)`` — stacking it on a foundation test is silent level
conflation.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from orpheus.diffusion.boundary_realizer import DiffusionBoundaryRealizer
from orpheus.diffusion.method_space import DiffusionMethodSpace
from orpheus.geometry.boundary import (
    AlbedoBoundary,
    BoundaryError,
    BoundaryTraceLaw,
    ConstantInflowSource,
    IdentityMap,
    IsotropicReturn,
    LambertianReemission,
    law_permutes_ordinates,
    PeriodicBoundary,
    PrescribedInflow,
    ReflectiveBoundary,
    ScalarResponse,
    SpecularReemission,
    SpecularReturn,
    VacuumInflow,
    WhiteBoundary,
    ZeroFluxBoundary,
)
from orpheus.numerics.operator import (
    IdentityOperator,
    ScaledOperator,
    ZeroOperator,
    adjointable,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.spaces.angular_trace_space import build_omega_dot_n
from orpheus.sn.boundary.realizer import SNBoundaryRealizer
from tests.sn._test_helpers import face_method_space


pytestmark = [pytest.mark.foundation]


# ═══════════════════════════════════════════════════════════════════════
# Fixtures — the four-quadrature list both routes share (see the header)
# ═══════════════════════════════════════════════════════════════════════

_FOUR_FACES = ("xmin", "xmax", "ymin", "ymax")

#: ⚠ **DO NOT TRIM THIS LIST.** The two routes' canonical mutations are
#: discriminated by DISJOINT quadratures (MEASURED 2026-08-01, table in this
#: module's header): the specular remap mutation reds ONLY on
#: ``gauss_legendre`` + ``lebedev17``, the diffuse ``> 0.0`` classification
#: mutation reds ONLY on ``product24``. Every row is either a discriminator
#: for one route or a documented CONTROL for the other; a "redundant
#: quadrature" cleanup makes one of the two anchors designed-green.
#: ``test_the_local_remap_is_not_the_identity_where_it_matters`` asserts this
#: split so it reds if the classification ever drifts.
#:
#: ``(id, quadrature_factory, face, axis, outward_sign, faces)``.
#:
#: ``gauss_legendre(8)`` is on the ``xmin`` face with ``outward_sign=-1`` so an
#: ``outward_sign`` hard-coded to ``+1`` reds; ``level_symmetric(6)`` is on
#: ``ymax`` with ``axis="y"`` so an ``axis`` hard-coded to ``"x"`` reds.
_FIXTURES = [
    ("gl4_xmax", lambda: Quadrature.gauss_legendre(4), "xmax", "x", +1,
     ("xmin", "xmax")),
    ("gl8_xmin", lambda: Quadrature.gauss_legendre(8), "xmin", "x", -1,
     ("xmin", "xmax")),
    ("product24_xmax", lambda: Quadrature.product(n_mu=2, n_phi=4), "xmax",
     "x", +1, _FOUR_FACES),
    ("lebedev17_xmax", lambda: Quadrature.lebedev(17), "xmax", "x", +1,
     _FOUR_FACES),
    ("ls6_ymax", lambda: Quadrature.level_symmetric(6), "ymax", "y", +1,
     _FOUR_FACES),
]

#: The three amplitude regimes are three DIFFERENT production branches:
#: ``_narrowed_zero_operator`` (α=0), ``ScaledOperator`` (0<α<1) and the bare
#: tensor product (α=1). A gate that covers only the middle one leaves the two
#: fast paths untested — which is exactly how the pre-B3.4b albedo arm kept two
#: un-narrowed endomorphisms alive at α ∈ {0, 1}.
_ALPHAS = [0.0, 0.3, 1.0]


def _space(fixture):
    _id, quad_of, face, _axis, _sign, faces = fixture
    quad = quad_of()
    return quad, face_method_space(quad, face=face, faces=faces)


def _probe(space, seed: int = 3) -> np.ndarray:
    r"""A random :math:`\Gamma_+` probe with DISTINCT trailing axes.

    ``(|Γ₊|, 4, 2)``: non-flat (a flat probe makes the Lambertian average equal
    its own input and nulls the whole redistribution), and ``4 != 2`` so a
    group ↔ spatial axis transpose cannot hide.
    """
    rng = np.random.default_rng(seed)
    return rng.normal(size=(int(space.outflow_indices.size), 4, 2))


def _sorted_half_traces(space) -> tuple[np.ndarray, np.ndarray]:
    return (
        np.sort(np.asarray(space.inflow_indices, dtype=np.intp)),
        np.sort(np.asarray(space.outflow_indices, dtype=np.intp)),
    )


# ═══════════════════════════════════════════════════════════════════════
# 1. The ≡ theorems — necessary, NOT sufficient
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("fixture", _FIXTURES, ids=[f[0] for f in _FIXTURES])
@pytest.mark.parametrize("alpha", _ALPHAS)
class TestEquivalenceTheorems:
    r"""``albedo(α, closure) ≡`` its geometry-tier sibling, as realized operators.

    **Tolerance: ``np.array_equal``, and the choice is structural.** Both routes
    reach ONE construction body with identical arguments, so the two operators
    hold the same permutation / kernel data and their ``apply`` executes the
    identical reduction tree. The predicted FP drift is EXACTLY zero
    (``vv-principles`` §bit-identity criterion 3: the bound is
    ``reduction_depth × ULP`` and the reduction-depth *difference* is 0). A
    non-zero drift means the routes did NOT execute the same body — precisely the
    bug this gate exists to catch — so ``assert_allclose`` at any rtol would
    admit the failure mode it is supposed to detect.

    **What this class canNOT catch.** A bug inside the shared body moves both
    sides identically. MEASURED: writing the ``to_local`` remap as a bare
    ``arange`` leaves every row here GREEN while the operator is wrong on three
    of five quadratures. The catchers are
    :class:`TestSpecularAgainstAnIndependentExpression` and
    :class:`TestDiffuseAgainstAnIndependentExpression`.
    """

    def test_specular_closure_equals_reflective(self, fixture, alpha) -> None:
        """T1 — forward action, bit-identical."""
        _id, _q, _face, axis, _sign, _faces = fixture
        quad, space = _space(fixture)
        x = _probe(space)
        realizer = SNBoundaryRealizer()
        got = realizer.realize(
            AlbedoBoundary(alpha, SpecularReturn(axis=axis)), space,
        ).apply(x)
        want = realizer.realize(
            ReflectiveBoundary(axis=axis, albedo=alpha), space,
        ).apply(x)
        np.testing.assert_array_equal(
            got, want,
            err_msg=(
                f"{_id} α={alpha}: albedo(α, SpecularReturn) and "
                f"reflective(axis, α) must realize to the SAME matrix — they "
                f"reach one construction body, so any difference is an "
                f"ARGUMENT drift (wrong axis / dropped α / wrong law_key)."
            ),
        )
        if alpha == 0.0:
            # α=0 is the zero map BY DESIGN; the non-vacuity leg below would
            # be false for the right reason, so assert the design instead.
            assert not np.any(want), "α=0 must be the zero map"
        else:
            assert np.count_nonzero(want), (
                f"{_id} α={alpha}: the comparison is 0 == 0 and proves nothing."
            )

    def test_isotropic_closure_equals_white(self, fixture, alpha) -> None:
        """T2 — forward action, bit-identical."""
        _id, _q, _face, axis, sign, _faces = fixture
        quad, space = _space(fixture)
        x = _probe(space)
        realizer = SNBoundaryRealizer()
        got = realizer.realize(
            AlbedoBoundary(alpha, IsotropicReturn(axis=axis, outward_sign=sign)),
            space,
        ).apply(x)
        want = realizer.realize(
            WhiteBoundary(axis=axis, outward_sign=sign, albedo=alpha), space,
        ).apply(x)
        np.testing.assert_array_equal(got, want, err_msg=f"{_id} α={alpha}")
        if alpha == 0.0:
            assert not np.any(want), "α=0 must be the zero map"
        else:
            assert np.count_nonzero(want), f"{_id} α={alpha}: vacuous"

    @pytest.mark.parametrize("closure_kind", ["specular", "isotropic"])
    def test_the_two_routes_agree_on_adjointability_and_on_the_transpose(
        self, fixture, alpha, closure_kind,
    ) -> None:
        r"""The capability set and the transpose action agree, both routes.

        The capability leg is separate from the value leg on purpose: a dropped
        ``CAP_APPLY_TRANSPOSE`` is SILENT at the value level (nothing raises
        until a consumer asks), and it gates the whole composite through
        :attr:`SNBoundaryOperator.is_adjointable`. So the claim is
        ``adjointable(a) == adjointable(b)`` — an EQUALITY, not
        ``is True``: MEASURED, the diffuse route is honestly NOT adjointable at
        α ∉ {0} (``LambertianReemission.is_adjointable`` is ``False`` until B5
        types it as the rank-one it is), and asserting ``True`` there would pin
        a fiction.
        """
        _id, _q, _face, axis, sign, _faces = fixture
        quad, space = _space(fixture)
        realizer = SNBoundaryRealizer()
        if closure_kind == "specular":
            a_law = AlbedoBoundary(alpha, SpecularReturn(axis=axis))
            b_law = ReflectiveBoundary(axis=axis, albedo=alpha)
        else:
            a_law = AlbedoBoundary(
                alpha, IsotropicReturn(axis=axis, outward_sign=sign),
            )
            b_law = WhiteBoundary(axis=axis, outward_sign=sign, albedo=alpha)
        a_op = realizer.realize(a_law, space)
        b_op = realizer.realize(b_law, space)
        assert adjointable(a_op) == adjointable(b_op), (
            f"{_id} α={alpha} {closure_kind}: the two routes disagree on "
            f"adjointability ({adjointable(a_op)} vs {adjointable(b_op)}) — a "
            f"capability drop is silent at the value level and gates the whole "
            f"SNBoundaryOperator composite."
        )
        if not adjointable(a_op):
            return
        rng = np.random.default_rng(7)
        y = rng.normal(size=(int(space.inflow_indices.size), 4, 2))
        got, want = a_op.apply_transpose(y), b_op.apply_transpose(y)
        assert got.shape == (int(space.outflow_indices.size), 4, 2), (
            f"{_id}: the transpose lands on Γ₊, so it must emit "
            f"{space.outflow_indices.size} rows; got {got.shape[0]}."
        )
        np.testing.assert_array_equal(got, want, err_msg=f"{_id} α={alpha}ᵀ")


# ═══════════════════════════════════════════════════════════════════════
# 2. The independent-expression anchors — the actual catchers
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("fixture", _FIXTURES, ids=[f[0] for f in _FIXTURES])
@pytest.mark.parametrize("alpha", _ALPHAS)
class TestSpecularAgainstAnIndependentExpression:
    r"""``albedo(α, SpecularReturn(a))`` vs a hand-written specular gather.

    The reference is four lines of numpy over ``quadrature.reflection_index``
    and the raw index sets — it shares nothing with the realizer above the
    ``np.take``. Row :math:`j` of the image reads the MIRROR of the
    :math:`j`-th inflow ordinate, which is an outflow ordinate, at its position
    inside :math:`\Gamma_+`.

    **Fixture discipline (MEASURED).** ``gauss_legendre(4/8)`` and
    ``lebedev(17)`` DISCRIMINATE the canonical remap mutation
    (``to_local`` → ``arange``); ``product(2,4)`` and ``level_symmetric(6)`` are
    BLIND to it, because there ``perm[sorted(inflow)] == sorted(outflow)`` and
    the local permutation IS ``arange``. The blind rows stay in the
    parametrisation as documented CONTROLS — they pin the coincidence, and a
    mutation run that reds the first three while leaving these two green is the
    proof the fixture list is load-bearing rather than decorative.

    ``np.array_equal``: a scale-then-gather introduces no re-association, so the
    predicted drift is 0 ULP.
    """

    def test_matches_the_hand_written_mirror_gather(self, fixture, alpha) -> None:
        _id, _q, _face, axis, _sign, _faces = fixture
        quad, space = _space(fixture)
        inflow, outflow = _sorted_half_traces(space)
        x = _probe(space)

        perm = quad.reflection_index(axis)
        full = np.zeros((int(quad.N),) + x.shape[1:], dtype=x.dtype)
        full[outflow] = x
        expected = alpha * full[perm[inflow]]

        got = SNBoundaryRealizer().realize(
            AlbedoBoundary(alpha, SpecularReturn(axis=axis)), space,
        ).apply(x)
        np.testing.assert_array_equal(
            got, expected,
            err_msg=(
                f"{_id} α={alpha}: the realized specular closure disagrees "
                f"with the hand-written mirror gather. This reference imports "
                f"nothing from the realizer, so a mismatch is a real defect in "
                f"the shared _specular_kernel body (which the ≡ theorem cannot "
                f"see)."
            ),
        )
        if alpha != 0.0:
            assert np.count_nonzero(expected), f"{_id}: vacuous comparison"

    def test_the_local_remap_is_not_the_identity_where_it_matters(
        self, fixture, alpha,
    ) -> None:
        r"""The fixture-discipline claim itself, asserted rather than trusted.

        Records — as a test, not a comment — WHICH quadratures discriminate the
        ``arange`` remap mutation. If a future quadrature change makes
        ``gauss_legendre`` agree with ``arange``, this reds and tells the reader
        the anchor above has gone blind, instead of letting it silently become
        theatre.
        """
        del alpha  # the remap is amplitude-independent
        _id, _q, _face, axis, _sign, _faces = fixture
        quad, space = _space(fixture)
        inflow, outflow = _sorted_half_traces(space)
        perm = quad.reflection_index(axis)
        local = np.searchsorted(outflow, perm[inflow])   # a LINEAR-scan-free
        naive = np.arange(inflow.size)                   # independent to_local
        discriminating = {"gl4_xmax", "gl8_xmin", "lebedev17_xmax"}
        if _id in discriminating:
            assert not np.array_equal(local, naive), (
                f"{_id} was the DISCRIMINATING fixture for the arange-remap "
                f"mutation and no longer is — the anchor above is now blind "
                f"here. Re-measure the fixture table in this module's header."
            )
        else:
            assert np.array_equal(local, naive), (
                f"{_id} was a documented BLIND control and now discriminates — "
                f"harmless, but the header table is stale."
            )


@pytest.mark.parametrize("fixture", _FIXTURES, ids=[f[0] for f in _FIXTURES])
@pytest.mark.parametrize("alpha", _ALPHAS)
class TestDiffuseAgainstAnIndependentExpression:
    r"""``albedo(α, IsotropicReturn(a, s))`` vs a hand-written Lambertian average.

    Reference: contract the :math:`\Gamma_+` rows against
    :math:`w\,|\Omega\cdot\hat n|`, divide by that weight's own sum over
    :math:`\Gamma_+`, broadcast over :math:`\Gamma_-`, scale by α. That
    normalisation — total returned current = α × incident current, NOT the
    textbook :math:`1/\pi` — is the convention
    :class:`~orpheus.sn.boundary.angular.AngularAverageOperator` implements and
    :class:`~orpheus.geometry.boundary.WhiteBoundary` documents.

    **Tolerance ``rtol=1e-14``, and it is the BOUND rather than a rounded-up
    measurement.** Production computes ``(cos_w[:,None,…] * psi).sum(axis=0)``;
    the reference computes ``np.tensordot(cos_w, psi, axes=(0,0))``. Same value,
    different reduction ORDER, so the bound is ``|Γ₊| × ULP`` — on
    ``lebedev(17)`` that is ``49 × 2.2e-16 ≈ 1.1e-14``. MEASURED max relative
    error: 5.1e-16 (gl4), 7.1e-16 (product24), 2.4e-16 (ls6), 3.2e-15
    (lebedev17) — inside the bound, largest on the deepest reduction, which is
    the predicted shape.

    **Fixture discipline (MEASURED).** ``product(2,4)`` is the ONLY fixture that
    discriminates the B3.4a classification twin (``> 0.0`` instead of
    ``> TANGENTIAL_EPS`` — it mis-admits ordinates 1 and 5 there and nothing
    anywhere else). It is the exact complement of the specular route's
    discriminating set, so neither list may be pruned in favour of the other.
    """

    def test_matches_the_hand_written_cosine_weighted_average(
        self, fixture, alpha,
    ) -> None:
        _id, _q, face, axis, sign, _faces = fixture
        quad, space = _space(fixture)
        inflow, outflow = _sorted_half_traces(space)
        x = _probe(space)

        omega_dot_n = build_omega_dot_n(quad, (face,))[0]
        out_w = (np.asarray(quad.weights, dtype=float) * np.abs(omega_dot_n))[
            outflow
        ]
        current = np.tensordot(out_w, x, axes=(0, 0))
        expected = np.broadcast_to(
            alpha * current / out_w.sum(), (inflow.size,) + current.shape,
        )

        got = SNBoundaryRealizer().realize(
            AlbedoBoundary(alpha, IsotropicReturn(axis=axis, outward_sign=sign)),
            space,
        ).apply(x)
        np.testing.assert_allclose(
            got, expected, rtol=1e-14, atol=0.0,
            err_msg=(
                f"{_id} α={alpha}: the realized Lambertian disagrees with the "
                f"hand-written cosine-weighted average beyond the "
                f"|Γ₊|×ULP reduction-order bound."
            ),
        )
        if alpha != 0.0:
            assert np.count_nonzero(expected), f"{_id}: vacuous comparison"

    def test_the_image_is_isotropic_over_gamma_minus(self, fixture, alpha) -> None:
        """Every Γ₋ row carries the SAME value — the defining Lambertian law.

        An intrinsic-property leg (this project's standing rule that a
        math-bearing type ships a test of its defining law): a re-emission that
        is not constant over the incoming hemisphere is not Lambertian, however
        well it matches an average.
        """
        _id, _q, _face, axis, sign, _faces = fixture
        quad, space = _space(fixture)
        got = SNBoundaryRealizer().realize(
            AlbedoBoundary(alpha, IsotropicReturn(axis=axis, outward_sign=sign)),
            space,
        ).apply(_probe(space))
        np.testing.assert_array_equal(
            got, np.broadcast_to(got[0][None, ...], got.shape),
            err_msg=f"{_id} α={alpha}: the re-emitted flux is not isotropic",
        )


# ═══════════════════════════════════════════════════════════════════════
# 3. The exactly-one-of-G-R invariant
# ═══════════════════════════════════════════════════════════════════════


def _g_is_trivial(geometry_map) -> bool:
    """``G`` fixes no geometry — the deck group's identity element."""
    return isinstance(geometry_map, IdentityMap)


def _r_is_trivial(response_kernel) -> bool:
    r"""``R = I`` — the law asserts NO physics.

    ``ScalarResponse(0.0)`` is deliberately NOT trivial: returning nothing is a
    constitutive statement (it is vacuum's entire content), not an absence of
    one. Trivial means the identity kernel and nothing else.
    """
    return (
        isinstance(response_kernel, ScalarResponse)
        and response_kernel.alpha == 1.0
    )


_B5_XFAIL = (
    "B5 — ReflectiveBoundary(axis, α<1) carries TWO non-trivial factors "
    "(G = SpecularMirror AND R = α·I), which the §11.1 ruling forbids: a "
    "symmetry plane is a quotient of the domain and cannot absorb. That object "
    "IS AlbedoBoundary(α, SpecularReturn(axis)) wearing the geometry costume, "
    "and B3.4b built the honest spelling. It is unreachable from a tag "
    "(_law_from_tag hard-codes albedo=1.0), so nothing production-facing "
    "depends on it. B5 retires ReflectiveBoundary's `albedo` parameter — "
    "DELETE this marker there; the XPASS(strict) failure is the point."
)

#: ``(id, law, expect_xfail)`` — every registry law with representative
#: amplitudes and both closures. ``AlbedoBoundary(1.0)`` bare is deliberately
#: ABSENT and pinned by its own positive test below; the completeness leg proves
#: the union is the whole registry.
_INVARIANT_ROWS = [
    ("vacuum", VacuumInflow(), False),
    ("reflective_a1", ReflectiveBoundary(axis="x", albedo=1.0), False),
    ("reflective_a07", ReflectiveBoundary(axis="x", albedo=0.7), True),
    ("white_a1", WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0), False),
    ("white_a03", WhiteBoundary(axis="x", outward_sign=+1, albedo=0.3), False),
    ("albedo_bare_a0", AlbedoBoundary(albedo=0.0), False),
    ("albedo_bare_a05", AlbedoBoundary(albedo=0.5), False),
    ("albedo_specular_a1",
     AlbedoBoundary(1.0, SpecularReturn(axis="x")), False),
    ("albedo_specular_a05",
     AlbedoBoundary(0.5, SpecularReturn(axis="x")), False),
    ("albedo_isotropic_a05",
     AlbedoBoundary(0.5, IsotropicReturn(axis="x", outward_sign=+1)), False),
    ("periodic", PeriodicBoundary(), False),
    ("prescribed",
     PrescribedInflow(source=ConstantInflowSource(value=1.0)), False),
    ("zero_flux", ZeroFluxBoundary(), False),
]


class TestExactlyOneFactorIsNonTrivial:
    r"""§11.1's illegal-states-unrepresentable law, over the WHOLE registry.

    > **EXACTLY ONE of** :math:`G`, :math:`R` **is non-trivial.**

    Its contrapositive is the useful direction: a law that asserts any physics
    at all has :math:`G = \mathrm{id}`. Reflective and periodic put everything
    in :math:`G` (they are quotients of the domain, asserting no physics);
    vacuum, white and albedo put everything in :math:`R` (they are surfaces,
    fixing no geometry).

    Why the registry and not "the tag-reachable laws"
    -------------------------------------------------

    Scoping to tag-reachable laws would SILENTLY DROP the one row §12.4 wants
    flagged: ``_law_from_tag`` hard-codes ``albedo=1.0`` for reflective, so
    ``ReflectiveBoundary(α<1)`` is not tag-reachable. And "tag-reachable" is
    method-dependent — ``AlbedoBoundary`` is absent from ``SNMesh``'s registry
    and present in the diffusion one — so the invariant's SCOPE would depend on
    which registry you consulted. The invariant is a property of the law
    ALGEBRA, not of a method's admission list.
    """

    @pytest.mark.parametrize(
        "law_id,law,expect_xfail", _INVARIANT_ROWS,
        ids=[r[0] for r in _INVARIANT_ROWS],
    )
    def test_exactly_one_of_g_r_is_non_trivial(
        self, law_id, law, expect_xfail, request,
    ) -> None:
        if expect_xfail:
            request.node.add_marker(
                pytest.mark.xfail(strict=True, reason=_B5_XFAIL)
            )
        G, R = law.geometry_map, law.response_kernel
        n_non_trivial = int(not _g_is_trivial(G)) + int(not _r_is_trivial(R))
        assert n_non_trivial == 1, (
            f"{law_id}: {n_non_trivial} non-trivial factors "
            f"(G={G!r} trivial={_g_is_trivial(G)}, "
            f"R={R!r} trivial={_r_is_trivial(R)}). §11.1: exactly one of G, R "
            f"is non-trivial — a quotient of the domain asserts no physics "
            f"(R = I), a surface fixes no geometry (G = id)."
        )

    def test_bare_albedo_at_unit_amplitude_carries_NEITHER_factor(self) -> None:
        r"""``AlbedoBoundary(1.0)`` bare is the registry's ONE degenerate row.

        ``G = IdentityMap`` and ``R = ScalarResponse(1.0) = I``: **zero**
        non-trivial factors, so it sits outside the parametrised invariant
        above. This is NOT a deferred defect and NOT an oversight — it is the
        degenerate factorisation a SCALAR trace forces, where ``G`` is the
        identity by dimension and ``α = 1`` means "return everything". It is
        exactly the under-determined spelling an ANGULAR method refuses
        (:class:`TestBareAlbedoIsRefusedBySN`), which is why it is pinned here
        positively rather than carried as a permanent xfail: an xfail with no
        landing phase is a red that never flips.

        §12.4 of the design names ``ReflectiveBoundary(α<1)`` as the invariant's
        counter-example and does not name this one. It is the second.
        """
        law = AlbedoBoundary(albedo=1.0)
        assert _g_is_trivial(law.geometry_map)
        assert _r_is_trivial(law.response_kernel)

    def test_the_law_inventory_is_complete(self) -> None:
        """Every registered law kind appears above, or as a named exception.

        Without this, a law added to the registry tomorrow escapes the
        invariant silently — and the deliberate exclusion of
        ``AlbedoBoundary(1.0)`` would be indistinguishable from an oversight.
        """
        covered = {type(law).__name__ for _id, law, _x in _INVARIANT_ROWS}
        registered = {
            type(BoundaryTraceLaw.create(key)).__name__
            for key in BoundaryTraceLaw.registry
        }
        missing = registered - covered
        assert not missing, (
            f"law kind(s) {sorted(missing)} are registered but absent from the "
            f"§11.1 G/R invariant gate — add them to _INVARIANT_ROWS, or as a "
            f"named positive exception with its reasoning."
        )


class TestTheClosureTierIsAmplitudeFree:
    r"""``α`` has exactly ONE home — the law — and the closure is a pure shape.

    ``amplitude`` is a :class:`BoundaryResponseKernel` Protocol member the
    DIFFUSION realizer reads, so the kernel must carry α; if the closure carried
    it too there would be two sources of one number. The closure is therefore
    amplitude-free and :meth:`ReemissionClosure.kernel` is the α-instantiation
    morphism.
    """

    @pytest.mark.parametrize("alpha", [0.0, 0.3, 1.0])
    def test_specular_return_instantiates_at_the_laws_alpha(self, alpha) -> None:
        law = AlbedoBoundary(alpha, SpecularReturn(axis="y"))
        R = law.response_kernel
        assert isinstance(R, SpecularReemission)
        assert R.amplitude == alpha
        assert R.axis == "y"
        assert R.is_zero == (alpha == 0.0)
        assert R.is_adjointable is True
        assert not hasattr(law.reemission, "alpha"), (
            "the closure must be amplitude-FREE — α lives on the law"
        )

    @pytest.mark.parametrize("alpha", [0.0, 0.3, 1.0])
    def test_isotropic_return_instantiates_at_the_laws_alpha(self, alpha) -> None:
        law = AlbedoBoundary(
            alpha, IsotropicReturn(axis="y", outward_sign=-1),
        )
        R = law.response_kernel
        assert isinstance(R, LambertianReemission)
        assert R.amplitude == alpha
        assert (R.axis, R.outward_sign) == ("y", -1)
        assert not hasattr(law.reemission, "alpha")

    def test_geometry_map_is_the_identity_for_every_closure(self) -> None:
        r"""``G = IdentityMap`` UNCONDITIONALLY — the ruling's core claim.

        In particular a ``SpecularReturn`` closure does NOT put a mirror in
        ``G``. That is the whole content of the 2026-08-01 ruling: the specular
        pairing of a WALL is constitutive, so it lives in ``R``; only a
        symmetry-plane's pairing is geometry. Mutating
        ``AlbedoBoundary.geometry_map`` to return ``SpecularMirror`` for a
        specular closure is the canonical RED for this class AND for the G/R
        invariant above.
        """
        for closure in (
            None,
            SpecularReturn(axis="x"),
            IsotropicReturn(axis="z", outward_sign=-1),
        ):
            law = AlbedoBoundary(0.5, closure)
            assert isinstance(law.geometry_map, IdentityMap), (
                f"closure={closure!r}: G must stay the deck group's identity"
            )


# ═══════════════════════════════════════════════════════════════════════
# 4. The refusal — both directions
# ═══════════════════════════════════════════════════════════════════════


class TestBareAlbedoIsRefusedBySN:
    r"""SN refuses ``AlbedoBoundary(α, reemission=None)``, naming both completions.

    The message assertions are LOAD-BEARING, not decoration. A bare
    ``pytest.raises(BoundaryError)`` is satisfied by ``_outflow_restriction``'s
    missing-face guard, by ``assert_realizable``, or by any other refusal on the
    path — a false green of exactly the shape ``vv-principles`` Mode-8 warns
    about. The three substring legs pin THIS raise.

    Equally load-bearing: the raising call is the PRODUCTION entry point
    ``SNBoundaryRealizer().realize``, never a locally-built exception (Mode-8
    class 5), and the positive control below proves an arm that refused
    EVERYTHING would not pass.
    """

    @pytest.mark.parametrize("alpha", [0.0, 0.5, 1.0])
    def test_refused_at_every_amplitude(self, alpha) -> None:
        r"""α coverage is mandatory: the three amplitudes were three DIFFERENT
        pre-B3.4b branches (``ZeroOperator`` / ``IdentityOperator`` /
        ``ScaledOperator``), so a refusal added only to the general branch would
        leave the two fast paths realizing full-face endomorphisms — which is
        precisely how albedo stayed un-narrowed at α ∈ {0, 1} through B3.2.
        """
        quad = Quadrature.gauss_legendre(8)
        space = face_method_space(quad, face="xmax")
        with pytest.raises(BoundaryError) as exc:
            SNBoundaryRealizer().realize(AlbedoBoundary(alpha), space)
        assert exc.value.law == "albedo", (
            f"the refusal must be attributed to the albedo law; got "
            f"{exc.value.law!r} — a different guard fired and this leg would "
            f"otherwise credit it as the angular-resolution refusal."
        )
        message = str(exc.value)
        if alpha == 0.0:
            # α = 0 is refused on a DIFFERENT argument, so it makes a
            # different promise — see the dedicated test below.
            return
        for token in ("SpecularReturn", "IsotropicReturn"):
            assert token in message, (
                f"the refusal must NAME both completions so the caller can act "
                f"on it; {token!r} missing from: {message}"
            )
        if alpha != 0.0:
            # KEYED TO α deliberately. At α ≠ 0 the defect really is a pairing
            # by array position — nothing says which outgoing direction feeds
            # which incoming one. At α = 0 that reasoning does NOT apply:
            # nothing returns, so no pairing is needed and the law is complete
            # on an angular trace. The refusal at α=0 stands on a DIFFERENT
            # argument (it is a ``VacuumInflow`` twin, and realizability must
            # not be discontinuous at one exact float), which is why the
            # message assertion is scoped rather than blanket — asserting the
            # under-determination wording at α=0 would pin a false reason.
            assert "ARRAY POSITION" in message.upper(), (
                f"the refusal must name the DEFECT (a pairing by array "
                f"position, carrying no geometry), not merely decline. "
                f"Got: {message}"
            )

    def test_the_alpha_zero_refusal_names_the_vacuum_twin_not_underdetermination(
        self,
    ) -> None:
        r"""``AlbedoBoundary(0.0)`` is refused for a DIFFERENT reason.

        At α = 0 nothing returns, so no pairing is needed and
        ``R = 0`` IS complete on an angular trace — the under-determination
        argument does not reach this row. The refusal is still right, on two
        other grounds: realizability must not be discontinuous at one exact
        float (``ScalarResponse.is_zero`` is documented as an EXACT compare
        precisely because "a near-zero albedo is a weak reflector, not a
        vacuum"), and letting it realize keeps a second spelling of
        :class:`VacuumInflow` alive on the SN arm.

        So the message must NAME the twin rather than recite the α ≠ 0
        reason, and — the discriminating half — it must NOT claim
        under-determination, which would send a reader looking for a missing
        pairing that is not missing.
        """
        quad = Quadrature.gauss_legendre(8)
        space = face_method_space(quad, face="xmax")
        with pytest.raises(BoundaryError) as exc:
            SNBoundaryRealizer().realize(AlbedoBoundary(0.0), space)
        assert exc.value.law == "albedo"
        message = str(exc.value)
        assert "VacuumInflow" in message, (
            f"the α=0 refusal must name the twin it defers to; got: {message}"
        )
        assert "under-determined" not in message.lower(), (
            f"α=0 is NOT under-determined — nothing returns, so no pairing is "
            f"needed. Reciting the α≠0 reason here sends the reader after a "
            f"missing pairing that is not missing. Got: {message}"
        )
        # CONTROL: α ≠ 0 DOES cite under-determination, so the negative leg
        # above is discriminating rather than vacuously true of every message.
        with pytest.raises(BoundaryError) as nonzero:
            SNBoundaryRealizer().realize(AlbedoBoundary(0.5), space)
        assert "under-determined" in str(nonzero.value).lower()

    @pytest.mark.parametrize("alpha", [0.0, 0.5, 1.0])
    @pytest.mark.parametrize("closure_kind", ["specular", "isotropic"])
    def test_positive_control_a_completed_law_realizes(
        self, alpha, closure_kind,
    ) -> None:
        """The CONTROL. Without it, an albedo arm that raised for EVERY input
        would pass the refusal legs above."""
        quad = Quadrature.gauss_legendre(8)
        space = face_method_space(quad, face="xmax")
        closure = (
            SpecularReturn(axis="x") if closure_kind == "specular"
            else IsotropicReturn(axis="x", outward_sign=+1)
        )
        op = SNBoundaryRealizer().realize(AlbedoBoundary(alpha, closure), space)
        image = op.apply(np.ones((int(space.outflow_indices.size), 3)))
        assert image.shape == (int(space.inflow_indices.size), 3)

    def test_an_unknown_closure_names_the_available_ones(self) -> None:
        r"""``ReemissionClosure`` is a ``runtime_checkable`` Protocol, so a
        third-party shape satisfies it structurally and reaches the realizer.
        The realizer must then refuse with a DISTINCT message naming what it
        does support — not fall through to the "no isinstance dispatch case"
        raise, which would blame the LAW rather than the closure.
        """
        class _UnknownReturn:
            def kernel(self, alpha):        # noqa: D102 - a structural stand-in
                return ScalarResponse(alpha)

        quad = Quadrature.gauss_legendre(8)
        space = face_method_space(quad, face="xmax")
        with pytest.raises(BoundaryError) as exc:
            SNBoundaryRealizer().realize(
                AlbedoBoundary(0.5, _UnknownReturn()), space,   # type: ignore[arg-type]
            )
        assert exc.value.law == "albedo"
        message = str(exc.value)
        assert "SpecularReturn" in message and "IsotropicReturn" in message


class TestBothClosuresRefuseAMisdeclaredAxisTheSameWay:
    r"""A closure declares its own ``axis``; the method space independently
    names the face. When they disagree, BOTH arms must refuse in the same
    vocabulary — an attributed
    :class:`~orpheus.geometry.boundary.BoundaryError` carrying ``law="albedo"``.

    This is a **symmetry** claim about the two arms, not two independent
    refusal claims, and it is the reason it needs its own class. Until this
    carve the diffuse arm cross-checked the declared orientation against the
    installation face and raised an attributed error, while the specular arm
    fell through to ``TraceRestrictionOperator.to_local`` and surfaced a raw
    ``ValueError`` about out-of-range rows — an implementation detail of a
    helper the caller never invoked, with no ``law`` attribution and no
    statement of what was actually wrong.

    The asymmetry pre-dates the closure for
    ``ReflectiveBoundary(axis="y")`` installed on an x-face, but nothing could
    *declare* an axis capable of disagreeing with its installation face until
    the closure introduced one, so it had no reachable consumer.

    A caller cannot route on an exception type it cannot predict, so "both
    arms refuse" is worth less than "both arms refuse the same way".
    """

    @pytest.mark.parametrize(
        "closure_kind", ["specular", "isotropic"],
    )
    def test_a_y_axis_closure_on_an_x_face_raises_an_attributed_error(
        self, closure_kind: str,
    ) -> None:
        quad = Quadrature.level_symmetric(6)
        space = face_method_space(
            quad, face="xmax", faces=("xmin", "xmax", "ymin", "ymax"),
        )
        closure = (
            SpecularReturn(axis="y") if closure_kind == "specular"
            else IsotropicReturn(axis="y", outward_sign=+1)
        )
        # ``BoundaryError`` NARROWLY: a bare ``Exception`` here would score the
        # raw ValueError this test exists to forbid as a pass (vv Mode 8,
        # class 6).
        with pytest.raises(BoundaryError) as exc:
            SNBoundaryRealizer().realize(AlbedoBoundary(0.7, closure), space)
        assert exc.value.law == "albedo", (
            f"{closure_kind}: refused, but unattributed — a consumer routing "
            f"on law= cannot see this one."
        )
        # The message must name the mismatch, not an index range: a reader
        # given only row numbers cannot tell a misdeclared axis from a broken
        # quadrature table.
        assert "'y'" in str(exc.value)

    def test_positive_control_the_matching_axis_realizes(self) -> None:
        """Without this, an arm that refused every closure would pass above."""
        quad = Quadrature.level_symmetric(6)
        space = face_method_space(
            quad, face="xmax", faces=("xmin", "xmax", "ymin", "ymax"),
        )
        for closure in (
            SpecularReturn(axis="x"),
            IsotropicReturn(axis="x", outward_sign=+1),
        ):
            op = SNBoundaryRealizer().realize(
                AlbedoBoundary(0.7, closure), space,
            )
            image = op.apply(np.ones((int(space.outflow_indices.size), 2)))
            assert image.shape == (int(space.inflow_indices.size), 2)


class TestDiffusionIsUnmovedAndClosureBlind:
    r"""The diffusion arm realizes ``AlbedoBoundary`` exactly as before.

    On a SCALAR trace there is one boundary degree of freedom,
    :math:`J^- = \alpha J^+`, so the bare amplitude IS the complete law and
    ``G`` is forced to the identity by dimension. ``BC("albedo", albedo=…)``
    means what it always meant.

    The second test is the one the design does not ask for and most needs: the
    closure must be **INERT** here. §12.3 asserts that the scalar trace has
    nothing for a closure to fix; nothing pins it. A future "read
    ``response_kernel`` instead of ``amplitude``" tidy-up would silently change
    diffusion answers, and this is the only gate that would catch it.
    """

    _PROBE = np.array([1.0, 2.0])

    @pytest.mark.parametrize(
        "alpha,expected_type,expected",
        [
            (0.0, ZeroOperator, [0.0, 0.0]),
            (0.5, ScaledOperator, [0.5, 1.0]),
            (1.0, IdentityOperator, [1.0, 2.0]),
        ],
    )
    def test_the_bare_law_realizes_exactly_as_before(
        self, alpha, expected_type, expected,
    ) -> None:
        op = DiffusionBoundaryRealizer().realize(
            AlbedoBoundary(albedo=alpha), DiffusionMethodSpace(),
        )
        assert isinstance(op, expected_type)
        np.testing.assert_array_equal(op.apply(self._PROBE), expected)

    @pytest.mark.parametrize("alpha", [0.0, 0.5, 1.0])
    @pytest.mark.parametrize(
        "closure",
        [
            None,
            SpecularReturn(axis="x"),
            IsotropicReturn(axis="x", outward_sign=+1),
            IsotropicReturn(axis="y", outward_sign=-1),
        ],
        ids=["bare", "specular_x", "isotropic_x", "isotropic_y_min"],
    )
    def test_the_closure_is_inert_on_a_scalar_trace(self, alpha, closure) -> None:
        realizer, space = DiffusionBoundaryRealizer(), DiffusionMethodSpace()
        got = realizer.realize(
            AlbedoBoundary(alpha, closure), space,
        ).apply(self._PROBE)
        want = realizer.realize(
            AlbedoBoundary(albedo=alpha), space,
        ).apply(self._PROBE)
        np.testing.assert_array_equal(
            got, want,
            err_msg=(
                f"α={alpha} closure={closure!r}: the diffusion arm must be "
                f"BLIND to the re-emission closure — a scalar trace has one "
                f"degree of freedom and no angular distribution to fix. A "
                f"difference here means the arm started dispatching on the "
                f"response kernel's TYPE, which silently changes every "
                f"BC('albedo', …) answer."
            ),
        )


# ═══════════════════════════════════════════════════════════════════════
# 5. The α = 0 zero map — a capability this carve OPENED
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize(
    "law_id,build_law",
    [
        ("reflective", lambda: ReflectiveBoundary(axis="x", albedo=0.0)),
        ("white", lambda: WhiteBoundary(axis="x", outward_sign=+1, albedo=0.0)),
        ("albedo_specular",
         lambda: AlbedoBoundary(0.0, SpecularReturn(axis="x"))),
        ("albedo_isotropic",
         lambda: AlbedoBoundary(0.0, IsotropicReturn(axis="x", outward_sign=+1))),
    ],
)
class TestZeroAmplitudeRealizesTheNarrowedZeroMap:
    r"""``α = 0`` on ANY kernel route is the honest zero map :math:`\Gamma_+ \to \Gamma_-`.

    This is a NET-NEW capability B3.4b opened and the design text does not
    mention. Before the fold into ``_attenuated_kernel_operator``,
    ``ReflectiveBoundary(axis, 0.0)`` and ``WhiteBoundary(…, 0.0)`` were LEGAL
    laws (α=0 satisfies every invariant, sub-Markov included) that could not be
    realized at all: they died in the numerics layer on
    ``ScaledOperator``'s zero-scalar refusal. Folding four routes into one body
    is what made that visible, and one answer fixed all four.

    A perfectly absorbing surface IS a vacuum, and now says so with the same
    object — with a working transpose, which the ``ScaledOperator`` path never
    reached.
    """

    def test_realizes_to_the_zero_map_with_both_spaces_named(
        self, law_id, build_law,
    ) -> None:
        quad = Quadrature.gauss_legendre(4)
        space = face_method_space(quad, face="xmax")
        n_in = int(space.inflow_indices.size)
        n_out = int(space.outflow_indices.size)
        op = SNBoundaryRealizer().realize(build_law(), space)
        assert isinstance(op, ZeroOperator), (
            f"{law_id} at α=0 must be the ZeroOperator, not a scaled kernel"
        )
        image = op.apply(np.full((n_out, 4, 2), 3.0))
        assert image.shape == (n_in, 4, 2)
        np.testing.assert_array_equal(image, 0.0)

    def test_the_transpose_lands_on_gamma_plus_and_annihilates(
        self, law_id, build_law,
    ) -> None:
        r"""The zero map's DEFINING law, on a NON-ZERO probe.

        ``|Γ₊| == |Γ₋|`` on every quadrature in the tree, so a shape assertion
        alone cannot tell the two-space zero map from an endomorphic ``0.0 * x``
        echo (``vv-principles`` Mode 12). Feeding a non-zero Γ₋ probe and
        requiring a zero Γ₊ image is the value leg that leaves that functional.
        """
        quad = Quadrature.gauss_legendre(4)
        space = face_method_space(quad, face="xmax")
        op = SNBoundaryRealizer().realize(build_law(), space)
        assert adjointable(op), f"{law_id} at α=0 must expose a transpose"
        probe = np.arange(
            int(space.inflow_indices.size) * 4 * 2, dtype=float,
        ).reshape(int(space.inflow_indices.size), 4, 2) + 1.0
        image = op.apply_transpose(probe)
        assert image.shape == (int(space.outflow_indices.size), 4, 2)
        np.testing.assert_array_equal(image, 0.0)

    def test_negative_control_the_scaled_path_still_refuses_a_zero_scalar(
        self, law_id, build_law,
    ) -> None:
        """The α=0 answer comes from the REALIZER's routing, not from a change
        in the numerics layer.

        Without this leg, the tests above would also pass if someone had
        "fixed" ``ScaledOperator`` to accept 0.0 — a different, much wider
        change with its own blast radius.
        """
        del law_id, build_law
        with pytest.raises(ValueError, match="zero scalar"):
            ScaledOperator(0.0, IdentityOperator())


# ═══════════════════════════════════════════════════════════════════════
# 6. The ≡ theorem at the FACTOR-READING consumers
# ═══════════════════════════════════════════════════════════════════════


class TestEquivalenceHoldsAtTheFactorReadingConsumers:
    r"""``≡`` must hold where production reads the FACTORS, not only where it
    reads the realized output.

    This is the generalisation of this module's headline caveat, and the one
    place in it with real teeth. Everywhere else the two routes share a
    realization body, so :class:`TestEquivalenceTheorems` cannot see a bug in
    it. Here they share **nothing**: four production sites answer *"does this
    law's realized* :math:`R\,G` *permute the angular index?"* by inspecting
    the law's factors, and ``ReflectiveBoundary`` carries its pairing in
    :math:`G` while ``AlbedoBoundary(α, SpecularReturn(a))`` carries it in
    :math:`R`. Two DIFFERENT interrogations of two DIFFERENT tiers must return
    the same answer for two laws that realize to the same matrix.

    They did not, before
    :func:`~orpheus.geometry.boundary.law_permutes_ordinates` existed. Each
    site spelled ``law.geometry_map.permutes_ordinates`` inline — correct
    while every specular pairing lived in :math:`G`, and half-right the moment
    B3.4b put one in :math:`R`. The consequences were not cosmetic:

    * the sweep schedule's reflecting-face set decides the
      ``B_lower`` / ``B_upper`` split — an unlagged face that couples ordinate
      :math:`n` to its mirror partner is a wrong fixed point, not a slow one;
    * the DSA low-order admission guard;
    * ``_has_ruled_corner_action`` — whether the off-quadrature
      :math:`\mu = \pm 1` ray is expressible;
    * ``RadialCharacteristicBoundaryOperator._reflect_corner`` — the swap.

    So the gate is stated as the equivalence itself (``albedo+specular``
    answers what ``reflective`` answers) rather than as four hard-coded
    booleans: a table of expected values would drift with the design, while
    the equivalence IS the design.

    .. note::

       ``SNMesh.BOUNDARY_OPERATOR_REGISTRY`` does not admit ``albedo`` yet
       (#189) and all four sites read ``sn_mesh.bc[face].law``, so this is a
       LATENT correction — it fires the day the law is registered. That is
       exactly why it needs a gate now: an unreachable predicate that is
       already false is a landmine, and nothing downstream would red.
    """

    #: ``(id, albedo_route, geometry_tier_sibling)`` — the pairs whose
    #: realized matrices are equal (pinned by ``TestEquivalenceTheorems``).
    _PAIRS = [
        ("specular",
         AlbedoBoundary(0.7, SpecularReturn(axis="x")),
         ReflectiveBoundary(axis="x", albedo=0.7)),
        ("isotropic",
         AlbedoBoundary(0.7, IsotropicReturn(axis="x", outward_sign=+1)),
         WhiteBoundary(axis="x", outward_sign=+1, albedo=0.7)),
    ]

    @pytest.mark.parametrize(
        "pair_id,albedo_law,sibling", _PAIRS, ids=[r[0] for r in _PAIRS],
    )
    def test_permutes_ordinates_agrees_across_the_two_tiers(
        self, pair_id, albedo_law, sibling,
    ) -> None:
        """Two laws that realize to the same matrix answer the same."""
        assert law_permutes_ordinates(albedo_law) == law_permutes_ordinates(
            sibling
        ), (
            f"{pair_id}: the albedo route and its geometry-tier sibling "
            f"realize to the SAME matrix but disagree on whether they permute "
            f"the angular index. Four production sites branch on this answer "
            f"— the sweep schedule's B_lower/B_upper split among them — so a "
            f"disagreement is two identical operators behaving differently."
        )

    def test_the_predicate_is_NOT_constant(self) -> None:
        """NON-VACUITY, and it is load-bearing here.

        The equality above is satisfied by a predicate that returns ``True``
        for everything, or ``False`` for everything — and ``False``-for-
        everything is precisely the pre-B3.4b bug's shape as seen from the
        albedo side. So pin both answers: the specular pair permutes, the
        diffuse pair does not.
        """
        assert law_permutes_ordinates(
            AlbedoBoundary(0.7, SpecularReturn(axis="x"))
        ) is True
        assert law_permutes_ordinates(
            AlbedoBoundary(0.7, IsotropicReturn(axis="x", outward_sign=+1))
        ) is False

    def test_the_bare_and_non_pairing_laws_do_not_permute(self) -> None:
        """SCOPE control: the predicate reads the PAIRING, not the tier.

        A law with no pairing in either tier must answer ``False`` whatever
        its amplitude — otherwise ``law_permutes_ordinates`` would be
        answering "is this an albedo law?" rather than the question it names,
        and the sweep schedule would lag faces that need no lagging.
        """
        for law in (
            VacuumInflow(),
            AlbedoBoundary(0.7),                      # bare — no closure
            WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0),
            PeriodicBoundary(),                       # spatial wrap, not angular
        ):
            assert law_permutes_ordinates(law) is False, (
                f"{law!r} carries no ANGULAR pairing in either tier"
            )

    def test_the_geometry_tier_alone_is_NOT_the_answer(self) -> None:
        r"""The retired inline spelling, pinned as the thing that must fail.

        ``law.geometry_map.permutes_ordinates`` — what all four sites used to
        read — returns ``False`` for ``AlbedoBoundary(α, SpecularReturn(a))``
        because that law's :math:`G` is honestly the identity. This asserts the
        DIVERGENCE, so the gate reds if someone "simplifies"
        :func:`law_permutes_ordinates` back to a single-tier read.
        """
        specular_wall = AlbedoBoundary(0.7, SpecularReturn(axis="x"))
        assert specular_wall.geometry_map.permutes_ordinates is False
        assert law_permutes_ordinates(specular_wall) is True, (
            "the whole point of the function: the pairing is in R here, so a "
            "geometry-tier-only read is exactly the half-right answer B3.4b "
            "created."
        )
