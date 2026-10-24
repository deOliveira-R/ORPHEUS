r"""Campaign phase **B3.4b** — the re-emission closure, and the discipline it enforces.

The user's 2026-08-01 ruling (`.claude/plans/archive/b3_domain_narrowing_crosswalk.md`
§11.1) put the specular pairing of an absorbing *surface* in :math:`R`, not
:math:`G`: ``AlbedoBoundary.geometry_map`` is ``SelfPairedDeck.identity()``
**unconditionally**, and a :class:`ReemissionClosure` selects the response
kernel. This module carries the four claims that ruling creates and that have no
older home:

1. **The ≡ theorems** (:class:`TestEquivalenceTheorems`) —
   ``albedo(α, SpecularReturn(a)) = α · reflective(a)`` (the wall is α times
   the mirror, which has no amplitude) and
   ``albedo(α, IsotropicReturn(a, s)) ≡ white(a, s, α)`` as realized operators.
2. **The independent-expression anchors**
   (:class:`TestSpecularAgainstAnIndependentExpression`,
   :class:`TestDiffuseAgainstAnIndependentExpression`) — the gates that actually
   catch a bug, because (1) cannot.
3. **The at-most-one-of-G-R invariant** (:class:`TestAtMostOneFactorIsNonTrivial`).
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
the independent geometric partner reference
(:func:`~tests._harness.references.mirror_partner_indices`) /
``quadrature.weights`` × ``omega_dot_n`` and import nothing from
:mod:`orpheus.sn.boundary.realizer` above the ``np.take`` /
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

RETIRED LEGS (wave 2, slice S3; the demotion check of each). The consumers these legs drove, ``SNBoundaryRealizer`` (S2b) and
``DiffusionBoundaryRealizer`` (S3), are retired, and so is the closure argument (the wall is ``WallCondition.specular(a)`` or
``.white(a)``, one type, no ``ReemissionClosure``). What each of the 16 legs claimed lives here now:

* the ≡ theorems (``TestEquivalenceTheorems``, 3 legs): the two spellings are ONE value (``specular(0)``, ``white(0)`` and
  ``vacuum()`` are equal, one digest: ``tests/gates/transport/test_wall_condition.py``), so the theorem is an identity of the type;
* the independent expressions (``TestSpecularAgainstAnIndependentExpression``, ``TestDiffuseAgainstAnIndependentExpression``, 3 legs):
  the hand-computed mirror gather, unit and partial, bit-exact, is ``tests/gates/sn/operators/test_sn_boundary_realizer.py``; the
  cosine-weighted average and the isotropy of the image are ``tests/gates/sn/operators/test_lambertian_chain.py``;
* the bare-albedo refusals and the misdeclared axis (``TestBareAlbedoIsRefusedBySN``, ``TestBothClosuresRefuseAMisdeclaredAxis...``,
  6 legs): unspellable, the wall has neither a bare albedo nor an axis;
* diffusion unmoved and closure-blind (``TestDiffusionIsUnmovedAndClosureBlind``, 2 legs): ``test_diffusion_vacuum_is_tier_1``
  and the albedo table of ``tests/gates/diffusion/test_boundary_realizer.py`` (a closure is not an argument of the wall);
* the zero amplitude (``TestZeroAmplitudeRealizesTheNarrowedZeroMap``, 2 legs): ``ZeroMap`` of
  ``tests/gates/sn/architecture/test_realized_boundary.py``.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from orpheus.geometry.boundary import (
    AlbedoBoundary,
    BoundaryTraceLaw,
    ConstantInflowSource,
    SelfPairedDeck,
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
from orpheus.numerics.operator import IdentityOperator, ZeroOperator
from orpheus.numerics.quadrature import Quadrature
from tests._harness.references import mirror_partner_indices
from tests.gates.sn._test_helpers import face_method_space


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

#: The three amplitude regimes were three DIFFERENT production branches until
#: S2 (``_narrowed_zero_operator`` at α=0, ``ScaledOperator`` for 0<α<1, the
#: bare tensor product at α=1); since S2 they are ONE spelling, ``α·(K_ω ⊗ I)``,
#: whose α=0 and α=1 members the algebra normalises (to the bound
#: ``ZeroMorphism`` and to ``K_ω ⊗ I`` itself). The three amplitudes stay the
#: rows: the two ends are where a normalisation defect would show — which is
#: exactly how the pre-B3.4b albedo arm kept two un-narrowed endomorphisms
#: alive at α ∈ {0, 1}.
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


# ═══════════════════════════════════════════════════════════════════════
# 2. The independent-expression anchors — the actual catchers
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("fixture", _FIXTURES, ids=[f[0] for f in _FIXTURES])
@pytest.mark.parametrize("alpha", _ALPHAS)
class TestSpecularAgainstAnIndependentExpression:
    r"""``albedo(α, SpecularReturn(a))`` vs a hand-written specular gather.

    The reference is four lines of numpy over the independent geometric
    partner map (``mirror_partner_indices``) and the raw index sets — it
    shares nothing with the realizer above the ``np.take``, and since §7d.3
    nothing with production's pairing derivation either. Row :math:`j` of the image reads the MIRROR of the
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
        perm = mirror_partner_indices(quad, axis)
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


# ═══════════════════════════════════════════════════════════════════════
# 3. The exactly-one-of-G-R invariant
# ═══════════════════════════════════════════════════════════════════════


def _g_is_trivial(geometry_map) -> bool:
    """``G`` fixes no geometry — the deck group's identity element."""
    return geometry_map == SelfPairedDeck.identity()


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


#: ``(id, law)`` — every registry law with representative amplitudes and
#: both closures. ``AlbedoBoundary(1.0)`` bare is deliberately
#: ABSENT and pinned by its own positive test below; the completeness leg proves
#: the union is the whole registry.
_INVARIANT_ROWS = [
    ("vacuum", VacuumInflow()),
    ("reflective_a1", ReflectiveBoundary(axis="x")),
    ("white_a1", WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)),
    ("white_a03", WhiteBoundary(axis="x", outward_sign=+1, albedo=0.3)),
    ("albedo_bare_a0", AlbedoBoundary(albedo=0.0)),
    ("albedo_bare_a05", AlbedoBoundary(albedo=0.5)),
    ("albedo_specular_a1",
     AlbedoBoundary(1.0, SpecularReturn(axis="x"))),
    ("albedo_specular_a05",
     AlbedoBoundary(0.5, SpecularReturn(axis="x"))),
    ("albedo_isotropic_a05",
     AlbedoBoundary(0.5, IsotropicReturn(axis="x", outward_sign=+1))),
    ("periodic", PeriodicBoundary()),
    ("prescribed",
     PrescribedInflow(source=ConstantInflowSource(value=1.0))),
    ("zero_flux", ZeroFluxBoundary()),
]


class TestAtMostOneFactorIsNonTrivial:
    r"""§11.1's illegal-states-unrepresentable law, over the WHOLE registry.

    > **AT MOST ONE of** :math:`G`, :math:`R` **is non-trivial.**

    (It read "exactly one" until 2026-10-01.) Every parametrized row below has
    EXACTLY one non-trivial factor; the one law with neither, the bare
    ``AlbedoBoundary(1.0)``, is its own test,
    :meth:`test_bare_albedo_at_unit_amplitude_carries_NEITHER_factor`.

    Its contrapositive is the useful direction: a law that asserts any physics
    at all has :math:`G = \mathrm{id}`. Reflective and periodic put everything
    in :math:`G` (they are quotients of the domain, asserting no physics);
    vacuum, white and albedo put everything in :math:`R` (they are surfaces,
    fixing no geometry).

    Why the registry and not "the tag-reachable laws"
    -------------------------------------------------

    Scoping to tag-reachable laws would have SILENTLY DROPPED the row §12.4
    wanted flagged, the attenuated mirror ``ReflectiveBoundary(α<1)`` (G a
    mirror AND R = α·I), which no tag could reach. It was a strict xfail here
    until the reflective cleanup (2026-10-01) retired the mirror's amplitude,
    so the counter-example can no longer be constructed (the gate is
    ``tests/gates/geometry/test_reflective_is_a_mirror.py``). And
    "tag-reachable" is
    method-dependent — ``AlbedoBoundary`` is absent from ``SNProblem``'s registry
    and present in the diffusion one — so the invariant's SCOPE would depend on
    which registry you consulted. The invariant is a property of the law
    ALGEBRA, not of a method's admission list.
    """

    @pytest.mark.parametrize(
        "law_id,law", _INVARIANT_ROWS,
        ids=[r[0] for r in _INVARIANT_ROWS],
    )
    def test_exactly_one_of_g_r_is_non_trivial(
        self, law_id, law,
    ) -> None:
        G, R = law.geometry_map, law.response_kernel
        n_non_trivial = int(not _g_is_trivial(G)) + int(not _r_is_trivial(R))
        assert n_non_trivial == 1, (
            f"{law_id}: {n_non_trivial} non-trivial factors "
            f"(G={G!r} trivial={_g_is_trivial(G)}, "
            f"R={R!r} trivial={_r_is_trivial(R)}). §11.1: at most one of G, R "
            f"is non-trivial, and every row here has exactly one — a quotient "
            f"of the domain asserts no physics (R = I), a surface fixes no "
            f"geometry (G = id)."
        )

    def test_bare_albedo_at_unit_amplitude_carries_NEITHER_factor(self) -> None:
        r"""``AlbedoBoundary(1.0)`` bare is the registry's ONE degenerate row.

        ``G = SelfPairedDeck.identity()`` and ``R = ScalarResponse(1.0) = I``: **zero**
        non-trivial factors, so it sits outside the parametrised invariant
        above. This is NOT a deferred defect and NOT an oversight — it is the
        degenerate factorisation a SCALAR trace forces, where ``G`` is the
        identity by dimension and ``α = 1`` means "return everything". It is
        exactly the under-determined spelling an ANGULAR method refuses
        (:class:`TestBareAlbedoIsRefusedBySN`), which is why it is pinned here
        positively rather than carried as a permanent xfail: an xfail with no
        landing phase is a red that never flips.

        §12.4 of the design named ``ReflectiveBoundary(α<1)`` as the
        invariant's counter-example (unconstructible since the reflective
        cleanup) and did not name this one. It was the second.
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
        covered = {type(law).__name__ for _id, law in _INVARIANT_ROWS}
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
        r"""``G = SelfPairedDeck.identity()`` UNCONDITIONALLY — the ruling's core claim.

        In particular a ``SpecularReturn`` closure does NOT put a mirror in
        ``G``. That is the whole content of the 2026-08-01 ruling: the specular
        pairing of a WALL is constitutive, so it lives in ``R``; only a
        symmetry-plane's pairing is geometry. Mutating
        ``AlbedoBoundary.geometry_map`` to return
        ``SelfPairedDeck.mirror(axis=...)`` for a
        specular closure is the canonical RED for this class AND for the G/R
        invariant above.
        """
        for closure in (
            None,
            SpecularReturn(axis="x"),
            IsotropicReturn(axis="z", outward_sign=-1),
        ):
            law = AlbedoBoundary(0.5, closure)
            assert law.geometry_map == SelfPairedDeck.identity(), (
                f"closure={closure!r}: G must stay the deck group's identity"
            )


# ═══════════════════════════════════════════════════════════════════════
# 4. The refusal — both directions
# ═══════════════════════════════════════════════════════════════════════


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


# ═══════════════════════════════════════════════════════════════════════
# 5. The α = 0 zero map — a capability this carve OPENED
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize(
    "law_id,build_law",
    [
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
    is what made that visible, and one answer fixed all four. The attenuated
    mirror's row left with its amplitude in the reflective cleanup
    (2026-10-01): a symmetry has no amplitude, and the specular wall at zero
    is the ``albedo_specular`` row.

    A perfectly absorbing surface IS a vacuum, and now says so with the same
    object — with a working transpose, which the ``ScaledOperator`` path never
    reached.
    """


    def test_the_numerics_layer_owns_the_zero_of_a_scaled_operator(
        self, law_id, build_law,
    ) -> None:
        """INVERTED at S2, by design (this row's own former docstring said S2
        would). Until S2 the α=0 answer came from the REALIZER's routing
        (``_narrowed_zero_operator``) and the numerics layer refused
        ``ScaledOperator(0.0, ·)``. Since S2 the zero is the algebra's
        (design D6, Q8): ``0·A`` is the zero map, so the realizer spells the
        α=0 member as ``0·(K_ω ⊗ I)`` like every other amplitude, and the rows
        above now witness the ALGEBRA's zero on the narrowed spaces. This row
        pins the numerics half: a spaceless pointwise operand's zero is the
        pointwise ``ZeroOperator()``; a bound operand's is ``ZeroMorphism`` on
        its own spaces (``test_linear_combination.py`` c09, c38).
        """
        del law_id, build_law
        assert isinstance(0.0 * IdentityOperator(), ZeroOperator)


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

       ``SNProblem.BOUNDARY_OPERATOR_REGISTRY`` does not admit the ``albedo``
       TAG (#189), but a typed law declared on the mesh reaches all four
       sites (``problem.bc[face].law``) since ``985497b5`` (2026-08-05). The
       correction was therefore live, not latent. Two of the four sites read
       a partial reflector as a mirror until 2026-09-30, because a "yes, it
       permutes" is not "it is a mirror": DSA's admission and the corner swap
       (ERR-094, whose third site, the eigenvalue's leakage term, sits outside
       this list). That is the case this gate was written for: a predicate
       believed unreachable is a landmine, and nothing downstream reds.
    """

    #: ``(id, albedo_route, geometry_tier_sibling)`` — the pairs whose
    #: realized matrices are equal (pinned by ``TestEquivalenceTheorems``).
    #: The specular pair is the perfect wall against the mirror, at α = 1:
    #: the mirror has no amplitude, so that is the only α at which the two
    #: tiers realize to one matrix (until the reflective cleanup the pair
    #: was at 0.7 against the attenuated mirror ``ReflectiveBoundary(x,
    #: 0.7)``, which no longer exists).
    _PAIRS = [
        ("specular",
         AlbedoBoundary(1.0, SpecularReturn(axis="x")),
         ReflectiveBoundary(axis="x")),
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
