r"""One body for every deck transformation's trace arrow (step G6.3-7).

Sibling of ``test_specular_deck_chain.py`` (the mirror's binding) and
``test_lambertian_chain.py`` (the two-link response chain). Step 7 replaces
``_specular_kernel``'s ``quadrature.reflection_index(axis)`` table lookup AND
the periodic arm's unbound ``IdentityOperator() & IdentityOperator()`` with a
single ``_deck_kernel`` that derives the arrow from the deck element's
**rigid motion**.

These gates were designed against the step-7b **contract**, before the
implementation existed, and are pinned to it rather than to what landed —
``scratch/step7_deck_uplift_verification_plan.md`` §3.2.

The contract, restated so a reader need not fetch the plan
----------------------------------------------------------

``h`` is the deck element carrying the DOMAIN face onto the INSTALLATION face,
:math:`h(f') = f`. Invariance of :math:`\psi` under the deck group gives

.. math::

   \gamma_-\psi|_f(\Omega) \;=\; \gamma_+\psi|_{f'}(h^{-1}\Omega),

so with :meth:`~orpheus.geometry.transformation.RigidMotion.permutes`'
convention :math:`h(x_i) = x_{\pi(i)}` the row for inflow ordinate ``j``
gathers global ordinate :math:`\pi^{-1}(j)`, remapped into
:math:`\Gamma_+(f')`'s local index by ``to_local``.

⭐⭐ Three measured facts that shaped these gates
------------------------------------------------

**1. The permutation is the LINEAR PART's, not the motion's.** `[M]`
2026-08-07: :meth:`RigidMotion.permutes` applies the **affine** action
(``on_points``), so for any deck element carrying a translation it moves the
nodes off the sphere and returns ``None``:

===========================  ==============  ==============
motion                       affine          linear part
===========================  ==============  ==============
mirror-x                     ``perm``        ``perm`` (equal)
mirror-x seated at 2.5       ⛔ ``None``     ``perm``
**wrap** ``t = ê_x``         ⛔ ``None``     ``perm`` = ``arange``
wrap ``t = 17 ê_x``          ⛔ ``None``     ``perm`` = ``arange``
rotation 90° about z         ``perm``        ``perm`` (equal)
rotation 90° **seated**      ⛔ ``None``     ``perm``
glide                        ⛔ ``None``     ``perm``
===========================  ==============  ==============

Measured on ``gauss_legendre(8)`` and ``product(4,8)``. The affine spelling is
correct **only** for a motion whose translation is exactly zero — which is
every element :class:`SelfPairedDeck` admits (its guard forbids the affine
part) and **no** wrap. So a kernel spelled ``motion.permutes(nodes)`` is
correct for the specular arm and hard-fails on the one arm step 7 exists to
build. :class:`TestThePermutationIsTheLinearParts` is that claim's gate; it is
posed on the OBSERVABLE (the kernel is invariant under the translation), so it
holds whichever way the fix is spelled — an inline ``RigidMotion(m.linear)`` or
a new ``permutes_directions`` on the type.

**2. The wrong gather direction is REFUSED, not silently wrong.** `[M]` on all
six admissible rotation fixtures, :math:`\pi^{-1}[\text{inflow}]` lands inside
:math:`\Gamma_+(f')` and :math:`\pi[\text{inflow}]` lands outside it — so
``to_local`` raises. Structural: :math:`h(\Gamma_-(f)) = h^2(\Gamma_+(f'))`,
which equals :math:`\Gamma_+(f')` **iff h is an involution** — exactly the case
(mirror, translation) where the two conventions coincide anyway. Consequence
for credit: the π-mutation's red must NOT be attributed to the value row
(a mutation that reds by RAISING has attributed nothing). The value row names
the right answer; :meth:`TestTheGatherIsThePreimage.test_the_forward_gather_is_OUT_of_range`
is where the discrimination is pinned.

**3. ``product(4,4)`` is the DEGENERATE rotation fixture.** `[M]` its local
permutation is exactly ``arange``, like the wrap's — so a keystone written
there is green under ``local := np.arange(n)``. The keystone runs on
``product(4,8)`` (``[1,2,0,4,5,3,7,8,6,10,11,9]``); ``product(4,4)`` is kept
as an explicitly-labelled degeneracy control.

Tolerances, measured per law rather than inherited
--------------------------------------------------

`[M]` the witness residual :math:`\max|x_{g} - h^{-1}x_j|`: **exactly 0.0** for
the mirror (a signed permutation is exact in IEEE) and for the wrap (its
``on_directions`` is the identity); **≤ 1.11e-16** for the 90° rotation, whose
matrix carries ``cos(π/2)`` round-off. So the mirror and wrap rows assert
``array_equal`` — a tolerance there would admit the bug — and only the rotation
rows carry ``atol``, set at ``1e-15`` from that measurement rather than at the
kernel's own ``1e-13`` match window.
"""

from __future__ import annotations

import dataclasses


import numpy as np
import pytest

from orpheus.geometry.boundary import (
    PairedDeck,
    PeriodicBoundary,
    ReflectiveBoundary,
    SelfPairedDeck,
)
from orpheus.geometry.boundary._errors import BoundaryError
from orpheus.geometry.transformation import RigidMotion
from orpheus.numerics.face_layout import FaceLayout
from orpheus.numerics.operator import adjointable
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.spaces import AngularTraceSpace
from orpheus.numerics.symmetry import _embedded_nodes
from orpheus.sn.boundary import realizer as realizer_module
from orpheus.sn.boundary.realizer import (
    SNBoundaryRealizer,
    _deck_kernel,
    _outflow_restriction,
)
from orpheus.sn.mesh.method_space import SNMethodSpace

pytestmark = pytest.mark.l1

_NG = 2

#: The measured witness tolerance for a rotation deck (see the module
#: docstring). The mirror and the wrap are asserted bit-exact and never use it.
_ROTATION_ATOL = 1e-15


# ---------------------------------------------------------------------------
# THE ONE CALL SITE
#
# Every gate below calls ``_kernel``; nothing calls ``_deck_kernel`` directly,
# so a signature change is a single edit here — the same Pattern-2 discipline
# the production carve applies to the two deck arms.
#
# ``crossed_diagnosis`` is the caller-supplied semantics for an ERR-045
# crossing. These gates never assert on its TEXT (each law's arm owns its own
# wording); they assert that the refusal arrives as a typed ``BoundaryError``.
# ---------------------------------------------------------------------------


def _kernel(quadrature, method_space, gamma_out_domain, *, motion,
            domain_face, law_key):
    """The production deck kernel, called with a gate-owned diagnosis."""
    return _deck_kernel(
        quadrature, method_space, gamma_out_domain,
        motion=motion,
        domain_face=domain_face,
        law_key=law_key,
        crossed_diagnosis=(
            f"gate fixture: {law_key} deck arrow "
            f"{domain_face} -> {method_space.face}"
        ),
    )


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

_X_FACES = ("xmin", "xmax")
_XY_FACES = ("xmin", "xmax", "ymin", "ymax")

#: ``name -> (factory, layout_faces)``.
#:
#: ``gauss_legendre(4)`` is deliberately absent: it carries no tangential
#: ordinates at all, so it nulls the very tier the half-traces separate. Same
#: rationale as ``test_specular_deck_chain.py``.
#:
#: The slab rule gets an x-only layout: its nodes are embedded as
#: :math:`(\mu, 0, 0)`, so every ordinate is TANGENTIAL to a y-face and the
#: half-traces there are empty — a layout carrying them would be modelling a
#: face the quadrature cannot resolve.
_QUADS: dict[str, tuple] = {
    "gauss_legendre(8)": (
        lambda: Quadrature.gauss_legendre(n_ordinates=8), _X_FACES,
    ),
    "product(2,4)": (lambda: Quadrature.product(n_mu=2, n_phi=4), _XY_FACES),
    "product(4,4)": (lambda: Quadrature.product(n_mu=4, n_phi=4), _XY_FACES),
    "product(4,5)": (lambda: Quadrature.product(n_mu=4, n_phi=5), _XY_FACES),
    "product(4,8)": (lambda: Quadrature.product(n_mu=4, n_phi=8), _XY_FACES),
    "product(4,12)": (lambda: Quadrature.product(n_mu=4, n_phi=12), _XY_FACES),
    "level_symmetric(6)": (
        lambda: Quadrature.level_symmetric(sn_order=6), _XY_FACES,
    ),
    "lebedev(17)": (lambda: Quadrature.lebedev(order=17), _XY_FACES),
}

#: Fixtures on which a specular mirror about x is realizable.
_MIRROR_CASES = [
    "gauss_legendre(8)", "product(2,4)", "product(4,4)", "product(4,8)",
    "level_symmetric(6)", "lebedev(17)",
]

#: Fixtures closed under a 90° rotation about z. `[M]` ``gauss_legendre(8)``
#: is NOT (its nodes lie on the x axis) and ``product(4,5)`` is NOT
#: (``n_phi % 4 != 0``); both serve as refusal witnesses instead.
_ROTATION_CASES = [
    "product(2,4)", "product(4,4)", "product(4,8)", "product(4,12)",
    "level_symmetric(6)", "lebedev(17)",
]

#: ⭐ Of those, the ones whose LOCAL permutation is not ``arange``. `[M]`
#: ``product(2,4)`` and ``product(4,4)`` are degenerate and are excluded here;
#: see ``test_the_small_product_rules_are_a_DEGENERACY_control``.
_ROTATION_NONTRIVIAL_LOCAL = [
    "product(4,8)", "product(4,12)", "level_symmetric(6)", "lebedev(17)",
]


def _trace(quad: Quadrature, faces) -> AngularTraceSpace:
    layout = FaceLayout.from_named_shapes(
        [(f, (int(quad.N), _NG)) for f in faces]
    )
    return AngularTraceSpace.from_quadrature_and_layout(quad, layout)


def _deck_fixture(case: str, *, face: str, domain_face: str):
    """``(quad, trace, method_space(face), gamma_out(domain_face), nodes)``.

    ``gamma_out`` is built on the DOMAIN face — which is the whole content of
    a genuine pair, and the thing a same-face restriction would silently get
    right for the mirror and wrong for everything else.
    """
    factory, faces = _QUADS[case]
    quad = factory()
    trace = _trace(quad, faces)
    method_space = SNMethodSpace.for_face(
        quadrature=quad, face=face, trace=trace
    )
    domain_space = SNMethodSpace.for_face(
        quadrature=quad, face=domain_face, trace=trace
    )
    gamma_out = _outflow_restriction(domain_space, "deck")
    return quad, trace, method_space, gamma_out, _embedded_nodes(quad.measure)


def _mirror(axis: str = "x") -> RigidMotion:
    return RigidMotion.reflection(normal=np.eye(3)[("x", "y", "z").index(axis)])


def _rotation_90() -> RigidMotion:
    r"""``R_z(+π/2)``: ``(x,y,z) ↦ (−y, x, z)``, so it carries ``ymin`` onto
    ``xmin`` — i.e. :math:`h(f') = f` with ``f' = ymin``, ``f = xmin``.

    Checked both ways in the plan (§5.1): the normals anti-correspond
    (:math:`h(\hat n_{f'}) = -\hat n_f`, the same relation the periodic pair
    satisfies, because a deck element maps the domain into the ADJACENT copy),
    and :math:`h^{-1}` carries every direction inflowing at ``xmin`` to one
    outflowing at ``ymin``.
    """
    return RigidMotion.rotation_about_axis(axis=[0.0, 0.0, 1.0],
                                           angle=np.pi / 2)


def _linear_part(motion: RigidMotion) -> RigidMotion:
    """The origin-fixing element with the same action on DIRECTIONS."""
    return RigidMotion(motion.linear)


def _global_gather(kernel, gamma_out) -> np.ndarray:
    """The kernel's LOCAL permutation lifted back to global ordinate indices.

    ``to_local`` is a bijection onto the restriction's index set, so this is
    exact and loses nothing — it just states the claim at the tier where an
    independent reference can be written.
    """
    return np.asarray(gamma_out.indices, dtype=np.intp)[np.asarray(kernel.perm)]


# ============================================================================
# I — the permutation is the LINEAR PART's  (the measured defect, §1 above)
# ============================================================================


class TestThePermutationIsTheLinearParts:
    r"""⭐⭐ A deck element's action on ORDINATES ignores its translation.

    :meth:`RigidMotion.on_points` and :meth:`RigidMotion.on_directions` are
    spelled separately precisely so that applying an affine map to a direction
    is unwriteable. :meth:`RigidMotion.permutes` is the one place that
    re-introduces it — it matches ``on_points(x)`` against the set — and the
    deck kernel is its first consumer whose motions carry translations.

    The gates are posed on the OBSERVABLE, not on the spelling: two motions
    with the same linear part must produce bit-identical kernels. That holds
    for an inline ``RigidMotion(m.linear)`` and for a future
    ``permutes_directions`` alike, and it reds for the affine spelling.
    """

    def test_the_AFFINE_spelling_refuses_the_wrap(self) -> None:
        r"""`[M]` The reason this class exists, asserted rather than narrated.

        A scope statement in the sibling module's idiom: if
        :meth:`RigidMotion.permutes` ever starts acting on directions, this
        row reds and the class's RATIONALE needs rewriting — not its
        assertions, which are about the kernel.
        """
        quad = Quadrature.gauss_legendre(n_ordinates=8)
        nodes = _embedded_nodes(quad.measure)
        wrap = PairedDeck.wrap(axis="x").motion
        assert wrap.permutes(nodes, atol=1e-13) is None, (
            "RigidMotion.permutes applies the AFFINE action, so a translation "
            "moves the ordinates off the sphere. If this now returns a "
            "permutation, `permutes` changed semantics and the deck kernel's "
            "linear-part extraction may have become redundant"
        )
        linear = _linear_part(wrap).permutes(nodes, atol=1e-13)
        assert linear is not None
        np.testing.assert_array_equal(linear.indices, np.arange(quad.N))

    @pytest.mark.parametrize("case", ["gauss_legendre(8)", "product(4,8)"])
    def test_a_SEATED_mirror_gives_the_same_kernel_as_an_origin_one(
        self, case: str
    ) -> None:
        r"""The mirror's offset is invisible — and here it must be *invisible*,
        not *unspellable*.

        :class:`SelfPairedDeck` closes the Mode-12 offset hole by REFUSING an
        affine element, so the mirror reaching this kernel is always
        origin-fixing. But the kernel is shared, and its contract is the
        general one: a seated element acts on directions exactly as its linear
        part does. Bit-exact — a signed permutation is exact in IEEE.
        """
        quad, _, method_space, gamma_out, _ = _deck_fixture(
            case, face="xmin", domain_face="xmin"
        )
        plain = _kernel(quad, method_space, gamma_out, motion=_mirror("x"),
                        domain_face="xmin", law_key="reflective")
        seated = _kernel(
            quad, method_space, gamma_out,
            motion=RigidMotion.reflection(normal=[1.0, 0.0, 0.0], offset=2.5),
            domain_face="xmin", law_key="reflective",
        )
        np.testing.assert_array_equal(plain.perm, seated.perm)

    @pytest.mark.parametrize("case", ["gauss_legendre(8)", "product(2,4)"])
    def test_the_wrap_LENGTH_is_invisible_to_the_kernel(
        self, case: str
    ) -> None:
        r"""⭐ ``t = ê_x`` and ``t = 17 ê_x`` give the same arrow, bit for bit.

        This is what licenses :meth:`PairedDeck.wrap` storing a UNIT
        translation: the physical extent is configuration, and the realization
        provably cannot see it. Under the affine spelling BOTH calls raise, so
        the row reds — attributed by the wrap's own error message.
        """
        quad, _, method_space, gamma_out, _ = _deck_fixture(
            case, face="xmin", domain_face="xmax"
        )
        unit = _kernel(
            quad, method_space, gamma_out,
            motion=PairedDeck.wrap(axis="x").motion,
            domain_face="xmax", law_key="periodic",
        )
        long = _kernel(
            quad, method_space, gamma_out,
            motion=RigidMotion.translation_by([17.0, 0.0, 0.0]),
            domain_face="xmax", law_key="periodic",
        )
        np.testing.assert_array_equal(unit.perm, long.perm)
        np.testing.assert_array_equal(
            unit.perm, np.arange(unit.n),
            err_msg=(
                "the wrap's linear part is the identity, so its arrow is the "
                "trivial relabelling between two DISTINCT index sets — the "
                "local permutation must be arange"
            ),
        )

    @pytest.mark.parametrize("case", _ROTATION_NONTRIVIAL_LOCAL)
    def test_a_SEATED_rotation_gives_the_same_kernel(self, case: str) -> None:
        r"""The latent #178 arm: a real sector deck is seated at the core axis.

        Its linear part is what acts on angle, so seating must be invisible —
        and on a fixture whose local permutation is NOT ``arange``, so the row
        compares a non-trivial array rather than two identities.
        """
        quad, _, method_space, gamma_out, _ = _deck_fixture(
            case, face="xmin", domain_face="ymin"
        )
        origin = _kernel(quad, method_space, gamma_out, motion=_rotation_90(),
                         domain_face="ymin", law_key="sector")
        seated = _kernel(
            quad, method_space, gamma_out,
            motion=_rotation_90().seated_at([1.0, 1.0, 0.0]),
            domain_face="ymin", law_key="sector",
        )
        np.testing.assert_array_equal(origin.perm, seated.perm)
        assert not np.array_equal(origin.perm, np.arange(origin.n)), (
            f"{case} must have a NON-trivial local permutation, or this row "
            f"compares two identities and proves nothing about seating"
        )


# ============================================================================
# A — (retired) the specular arm's inheritance oracle
#
# ``TestTheSpecularArmInheritsTheRetiredTable`` lived here until §7d.3: a
# one-time proof that the motion-derived arrow equalled the
# ``reflection_index``-derived one (`[M]` bit-identical on 24/24
# quadrature × axis pairs, including the shared product(4,5) omission).
# It died WITH its oracle, per its own docstring — its surviving successor
# is ``TestTheGatherIsThePreimage``'s witness-check below, which verifies
# the claimed partner's COORDINATES instead of comparing derivations.
# The binding claim its second row carried (Γ₊(f) → Γ₋(f), by ``is``) is
# carried forward split by leg: the DOMAIN leg by the self-paired/paired
# split gate at the end of this module (``mirrored.domain is
# trace.outflow_space("xmin")``); the CODOMAIN leg rides on
# ``TestPeriodicIsBoundToThePartner`` alone — the split gate's codomain
# assertion is ``mirrored.codomain is realized.codomain``, which a
# None/None regression satisfies (`[M]` qa M5: dropping ``codomain=`` in
# the kernel reds 5 periodic rows and leaves the split gate green).
# ============================================================================


# ============================================================================
# B / C — the gather is the motion's PREIMAGE
# ============================================================================


class TestTheGatherIsThePreimage:
    r"""⭐⭐ The direction of the arrow, and the module's keystone.

    Two independent references, deliberately:

    * the **witness-check** — ``nodes[gather] == h⁻¹·nodes[inflow]``. It
      *verifies* the claimed partner instead of *searching* for one, so it
      shares no machinery with :meth:`RigidMotion.permutes` (a
      nearest-neighbour search) and survives ``reflection_index``'s
      retirement. Works for every law, including the ones where the direction
      is unobservable.
    * the **φ-lattice reference** — for a product rule the expected preimage is
      pure arithmetic on ``(level, p)``, with no search, no set membership and
      no call into the code under test. The strongest class available for this
      object (exact integers, no reference needed).
    """

    # -- the witness-check, every law -------------------------------------

    @pytest.mark.parametrize("case", _MIRROR_CASES)
    def test_the_MIRROR_gather_carries_the_reflected_direction(
        self, case: str
    ) -> None:
        """`[M]` residual exactly 0.0 — a signed permutation is IEEE-exact."""
        quad, _, method_space, gamma_out, nodes = _deck_fixture(
            case, face="xmin", domain_face="xmin"
        )
        motion = _mirror("x")
        kernel = _kernel(quad, method_space, gamma_out, motion=motion,
                         domain_face="xmin", law_key="reflective")
        inflow = np.asarray(method_space.inflow_indices, dtype=np.intp)
        gather = _global_gather(kernel, gamma_out)
        np.testing.assert_array_equal(
            nodes[gather],
            _linear_part(motion).inverse().on_directions(nodes[inflow]),
        )

    @pytest.mark.parametrize(
        "case", ["gauss_legendre(8)", "product(2,4)", "product(4,8)",
                 "level_symmetric(6)", "lebedev(17)"]
    )
    def test_the_WRAP_gather_carries_the_same_direction(
        self, case: str
    ) -> None:
        r"""`[M]` residual exactly 0.0 — a wrap's ``on_directions`` IS the
        identity, so the partner ordinate must carry the identical vector.

        This is the quotient reading as an assertion: ordinate ``n`` at the
        partner face feeds ordinate ``n`` here, at the SAME angle.
        """
        quad, _, method_space, gamma_out, nodes = _deck_fixture(
            case, face="xmin", domain_face="xmax"
        )
        kernel = _kernel(
            quad, method_space, gamma_out,
            motion=PairedDeck.wrap(axis="x").motion,
            domain_face="xmax", law_key="periodic",
        )
        inflow = np.asarray(method_space.inflow_indices, dtype=np.intp)
        gather = _global_gather(kernel, gamma_out)
        np.testing.assert_array_equal(nodes[gather], nodes[inflow])

    @pytest.mark.parametrize("case", _ROTATION_CASES)
    def test_the_ROTATION_gather_carries_the_inverse_rotated_direction(
        self, case: str
    ) -> None:
        r"""⭐ ``nodes[g] == R_z(−π/2)·nodes[j]``, ``atol`` from measurement.

        `[M]` max residual ``1.11e-16`` over all six fixtures; the tolerance
        is ``1e-15``, from that measurement — NOT the kernel's ``1e-13`` match
        window, which is nine orders looser than the code achieves and would
        accept a genuinely mismatched node on a coarse rule.

        This is the row that catches a gather in the FORWARD direction *by
        value* on any fixture where the forward gather happened to stay in
        range; on the shipped fixtures the range guard fires first, which is
        why the discrimination has its own row below.
        """
        quad, _, method_space, gamma_out, nodes = _deck_fixture(
            case, face="xmin", domain_face="ymin"
        )
        motion = _rotation_90()
        kernel = _kernel(quad, method_space, gamma_out, motion=motion,
                         domain_face="ymin", law_key="sector")
        inflow = np.asarray(method_space.inflow_indices, dtype=np.intp)
        gather = _global_gather(kernel, gamma_out)
        np.testing.assert_allclose(
            nodes[gather],
            motion.inverse().on_directions(nodes[inflow]),
            rtol=0.0, atol=_ROTATION_ATOL,
        )

    # -- the φ-lattice reference, product rules only ----------------------

    @pytest.mark.parametrize("case", ["product(4,4)", "product(4,8)",
                                      "product(4,12)"])
    def test_the_lattice_CONVENTION_this_class_relies_on(
        self, case: str
    ) -> None:
        r"""`[M]` ``index = level·n_φ + p``, ``μ = leggauss(n_mu)[level]``,
        ``φ = 2πp/n_φ`` — verified to ``< 1e-14``.

        Without this row the symbolic reference below is an assumption wearing
        the authority of arithmetic. With it, a change to the product rule's
        node ordering reds HERE first and says which assumption moved, instead
        of reddening the keystone with a misleading diagnosis.
        """
        n_mu, n_phi = _parse_product(case)
        quad = _QUADS[case][0]()
        nodes = _embedded_nodes(quad.measure)
        mu_nodes = np.polynomial.legendre.leggauss(n_mu)[0]
        for level in range(n_mu):
            mu = mu_nodes[level]
            sin_theta = np.sqrt(1.0 - mu * mu)
            for p in range(n_phi):
                phi = 2.0 * np.pi * p / n_phi
                np.testing.assert_allclose(
                    nodes[level * n_phi + p],
                    [sin_theta * np.cos(phi), sin_theta * np.sin(phi), mu],
                    rtol=0.0, atol=1e-14,
                )

    @pytest.mark.parametrize("case", ["product(4,8)", "product(4,12)"])
    def test_the_gather_is_the_PHI_LATTICE_preimage(self, case: str) -> None:
        r"""⭐⭐ **The keystone.** ``h⁻¹`` is ``φ ↦ φ − π/2``, so the source of
        ``index(level, p)`` is ``index(level, (p − n_φ/4) mod n_φ)``.

        The reference is built from ``(n_mu, n_phi)`` by integer arithmetic —
        it never calls :meth:`RigidMotion.permutes`, never calls
        ``reflection_index``, never touches the kernel. Structural
        independence in the strongest form this object admits (``vv`` L11 /
        lessons §8: exact integer arithmetic needs no reference at all).

        Run on the two fixtures whose local permutation is NOT ``arange``;
        ``product(4,4)`` is the degeneracy control below.
        """
        n_mu, n_phi = _parse_product(case)
        assert n_phi % 4 == 0, (
            f"{case}: a quarter turn maps the azimuth grid onto itself only "
            f"when 4 | n_phi; this fixture cannot host the keystone"
        )
        quad, _, method_space, gamma_out, _ = _deck_fixture(
            case, face="xmin", domain_face="ymin"
        )
        kernel = _kernel(quad, method_space, gamma_out, motion=_rotation_90(),
                         domain_face="ymin", law_key="sector")
        inflow = np.asarray(method_space.inflow_indices, dtype=np.intp)
        level, p = np.divmod(inflow, n_phi)
        expected_global = level * n_phi + (p - n_phi // 4) % n_phi
        np.testing.assert_array_equal(
            _global_gather(kernel, gamma_out), expected_global,
            err_msg=(
                "the gather must be the h^-1 PREIMAGE (phi -> phi - pi/2). "
                "Gathering pi[inflow] instead of pi^-1[inflow] shifts the "
                "azimuth the wrong way; a level/azimuth index transposition "
                "or a sign error in the angle move it differently again."
            ),
        )

    @pytest.mark.parametrize("case", _ROTATION_NONTRIVIAL_LOCAL)
    def test_the_LOCAL_remap_is_not_arange(self, case: str) -> None:
        r"""The ``to_local`` step, pinned separately from the gather.

        A correct global gather with a dropped ``to_local`` — the
        ``arange(g.size)`` shortcut ``to_local``'s own docstring warns about —
        reds here and only here. The reference is ``searchsorted`` against the
        restriction's sorted index set, written out rather than delegated, so
        the two sides are different expressions.
        """
        quad, _, method_space, gamma_out, _ = _deck_fixture(
            case, face="xmin", domain_face="ymin"
        )
        kernel = _kernel(quad, method_space, gamma_out, motion=_rotation_90(),
                         domain_face="ymin", law_key="sector")
        gather = _global_gather(kernel, gamma_out)
        expected_local = np.searchsorted(
            np.sort(np.asarray(gamma_out.indices)), gather
        )
        np.testing.assert_array_equal(kernel.perm, expected_local)
        assert not np.array_equal(kernel.perm, np.arange(kernel.n)), (
            f"{case} was chosen because its local permutation is NOT the "
            f"identity; if it has become the identity the fixture no longer "
            f"discriminates a dropped to_local and this list must be re-picked"
        )

    @pytest.mark.parametrize("case", ["product(2,4)", "product(4,4)"])
    def test_the_small_product_rules_are_a_DEGENERACY_control(
        self, case: str
    ) -> None:
        r"""⛔ **This row proves nothing about the local remap — that is its
        point.**

        `[M]` on ``product(2,4)`` and ``product(4,4)`` the rotation deck's
        local permutation is exactly ``arange``: one azimuth per μ-level lands
        in each half-trace, and both index sets sort the same way. A keystone
        written on these fixtures — which is what the step's brief proposed —
        is green under ``local := np.arange(n)``.

        Committed so the degeneracy is a recorded fact rather than a trap: if
        it ever stops holding, the keystone's fixture rationale has moved and
        this row says so.
        """
        quad, _, method_space, gamma_out, _ = _deck_fixture(
            case, face="xmin", domain_face="ymin"
        )
        kernel = _kernel(quad, method_space, gamma_out, motion=_rotation_90(),
                         domain_face="ymin", law_key="sector")
        np.testing.assert_array_equal(kernel.perm, np.arange(kernel.n))

    # -- the discrimination itself ----------------------------------------

    @pytest.mark.parametrize("case", _ROTATION_CASES)
    def test_the_forward_gather_is_OUT_of_range(self, case: str) -> None:
        r"""⭐ Where the π-vs-π⁻¹ mutation's red is honestly attributed.

        `[M]` on all six rotation fixtures ``π[inflow]`` lands in
        :math:`\Gamma_-(f')`, not :math:`\Gamma_+(f')`, so the domain
        space's
        :meth:`~orpheus.numerics.spaces.angular_trace_space.AngularFaceTraceSpace.to_local`
        refuses it. Structural:
        :math:`h(\Gamma_-(f)) = h^2(\Gamma_+(f'))`, and that equals
        :math:`\Gamma_+(f')` **iff h² = e** — i.e. only for the involutions,
        where the two conventions coincide anyway.

        Two jobs. It is the **activation guard** for the keystone (if the two
        gathers ever coincided, the keystone would be vacuous and this row
        reds first), and it pins that the wrong convention is *unrealizable*
        rather than silently wrong — free correctness insurance for the latent
        sector consumer (#178).
        """
        quad, _, method_space, gamma_out, nodes = _deck_fixture(
            case, face="xmin", domain_face="ymin"
        )
        pi = _linear_part(_rotation_90()).permutes(nodes, atol=1e-13)
        assert pi is not None
        inflow = np.asarray(method_space.inflow_indices, dtype=np.intp)
        forward = pi.indices[inflow]
        backward = pi.inverse().indices[inflow]
        assert not np.array_equal(forward, backward), (
            f"{case}: the two gather conventions coincide here, so every "
            f"row in this class is blind to the direction. A quarter turn is "
            f"not an involution — if this fires, the fixture is not the "
            f"motion it claims to be"
        )
        with pytest.raises(ValueError, match="are not ordinates of"):
            gamma_out.codomain.to_local(forward)
        # …and the correct one is accepted, so the refusal above is about the
        # DIRECTION and not about the space refusing everything.
        gamma_out.codomain.to_local(backward)


def _parse_product(case: str) -> tuple[int, int]:
    """``"product(4,8)" -> (4, 8)``."""
    inner = case[case.index("(") + 1: case.index(")")]
    n_mu, n_phi = (int(t) for t in inner.split(","))
    return n_mu, n_phi


# ============================================================================
# D / H — the refusals
# ============================================================================


class TestTheDeckKernelRefuses:
    r"""A motion that does not permute the ordinates must raise, attributed.

    ⚠ The message this inherits is a real hazard. The retired path raised a
    bare ``ValueError`` ("no precomputed reflection partner"), and
    :class:`~orpheus.geometry.boundary._errors.BoundaryError` **subclasses**
    ``ValueError`` — so a ``pytest.raises(ValueError)`` gate stays green across
    the whole change and proves nothing. Every row here names
    ``BoundaryError`` explicitly and matches on the message.
    """

    @pytest.mark.catches("ERR-074")
    def test_an_UNCERTIFIED_mirror_is_refused(self) -> None:
        r"""``product(4, 5)``: odd :math:`n_\varphi`, so :math:`\sigma_x` is
        not a symmetry of the azimuth grid.

        ERR-074's own site: the pre-2026-08-02 table returned a bare
        nearest-neighbour ``argmin`` here, wrong by ``0.42`` in the direction
        cosines and *still an involution*, so a self-inverse check could not
        see it. `[M]` both the certified table and
        :meth:`RigidMotion.permutes` now decline this axis — the refusal
        migrates 1:1, and this row is what proves the boundary tier kept it.
        """
        quad, _, method_space, gamma_out, nodes = _deck_fixture(
            "product(4,5)", face="xmin", domain_face="xmin",
        )
        assert _mirror("x").permutes(nodes, atol=1e-13) is None, (
            "product(4,5) must NOT be sigma_x-closed, or this row is not the "
            "refusal witness it claims to be — the mirror planes of a product "
            "rule sit at k*pi/n_phi and sigma_x needs one at pi/2"
        )
        with pytest.raises(BoundaryError) as excinfo:
            _kernel(quad, method_space, gamma_out, motion=_mirror("x"),
                    domain_face="xmin", law_key="reflective")
        assert excinfo.value.law == "reflective"

    def test_a_motion_the_quadrature_is_not_CLOSED_under_is_refused(
        self,
    ) -> None:
        r"""A 90° rotation on a slab quadrature: its nodes lie on the x axis,
        so :math:`R_z` carries them off the set entirely.

        A second refusal witness with a *different reason* from the row above
        (a whole family missing, versus a parity mismatch in one grid), so a
        single over-broad ``match=`` cannot cover both.
        """
        quad, _, method_space, gamma_out, _ = _deck_fixture(
            "gauss_legendre(8)", face="xmin", domain_face="xmin"
        )
        with pytest.raises(BoundaryError):
            _kernel(quad, method_space, gamma_out, motion=_rotation_90(),
                    domain_face="xmin", law_key="sector")

    def test_the_POSITIVE_control_the_same_call_succeeds(self) -> None:
        """Without this, a kernel that refused everything would pass above."""
        quad, _, method_space, gamma_out, _ = _deck_fixture(
            "product(4,4)", face="xmin", domain_face="xmin"
        )
        kernel = _kernel(quad, method_space, gamma_out, motion=_mirror("x"),
                         domain_face="xmin", law_key="reflective")
        assert kernel.n == gamma_out.n_restricted

    @pytest.mark.catches("ERR-073")
    def test_a_DUPLICATED_node_makes_the_match_non_injective_and_is_refused(
        self,
    ) -> None:
        r"""ERR-073's minimal witness, at the boundary tier.

        Appending a bit-identical duplicate of a node and of its mirror
        partner leaves the set closed under :math:`\sigma_x` as a *multiset*
        and leaves the two half-traces the SAME size — `[M]` 5 and 5 on
        ``gauss_legendre(8)`` + 2 duplicates — so the kernel's extent guard
        cannot fire first and the refusal is attributable to bijectivity
        alone. That is the clause a nearest-neighbour loop omits: every image
        has *a* partner, two sources land on one target, and one node is left
        unmatched.

        No tolerance games: the duplicate is bit-identical, so a looser match
        window cannot rescue it.
        """
        base = Quadrature.gauss_legendre(n_ordinates=8)
        nodes = np.asarray(base.measure.nodes)
        weights = np.asarray(base.measure.weights)
        duplicated = dataclasses.replace(
            base.measure,
            nodes=np.concatenate([nodes, nodes[[0, 7]]]),
            weights=np.concatenate([weights, weights[[0, 7]]]),
        )
        quad = Quadrature(measure=duplicated)
        nodes3 = _embedded_nodes(quad.measure)
        assert _mirror("x").permutes(nodes3, atol=1e-13) is None, (
            "the duplicated set must defeat the bijectivity clause, or this "
            "row is not the ERR-073 witness it claims to be"
        )
        assert _linear_part(_mirror("x")).permutes(
            _embedded_nodes(base.measure), atol=1e-13
        ) is not None, "the UNduplicated control must still permute"

        trace = _trace(quad, _X_FACES)
        method_space = SNMethodSpace.for_face(
            quadrature=quad, face="xmin", trace=trace
        )
        gamma_out = _outflow_restriction(method_space, "reflective")
        inflow = np.asarray(method_space.inflow_indices, dtype=np.intp)
        assert gamma_out.n_restricted == inflow.size, (
            "the two half-traces must stay the same size, or the extent guard "
            "fires first and the refusal is attributed to the wrong clause"
        )
        with pytest.raises(BoundaryError):
            _kernel(quad, method_space, gamma_out, motion=_mirror("x"),
                    domain_face="xmin", law_key="reflective")


# ============================================================================
# E / F — periodic, through the PUBLIC realizer
# ============================================================================

_PERIODIC_CASES = ["gauss_legendre(8)", "product(2,4)", "level_symmetric(6)",
                   "lebedev(17)"]


def _realized_periodic(case: str, face: str = "xmin"):
    """``(trace, realized_operator, partner)`` via the public realizer."""
    factory, faces = _QUADS[case]
    quad = factory()
    trace = _trace(quad, faces)
    method_space = SNMethodSpace.for_face(
        quadrature=quad, face=face, trace=trace
    )
    law = PeriodicBoundary(axis="x")
    realized = SNBoundaryRealizer().realize(law, method_space)
    return trace, realized, law.geometry_map.domain_face(face)


class TestPeriodicIsBoundToThePartner:
    r"""⭐ Periodic's link becomes a typed arrow — the step's stated OUTCOME.

    Every other law's angular factor is already one (``PermutationOperator``,
    ``OperatorProduct``, ``ZeroOperator``); periodic alone realized to
    ``IdentityOperator``, an **endomorphism** :math:`V \to V` that names no
    spaces and therefore can never be the isomorphism
    :math:`\Gamma_+(f') \to \Gamma_-(f)` a torus wrap needs.

    ⭐ These rows are the **only possible catcher** for a ``domain``↔
    ``codomain`` swap. Measured at step 5 on the specular sibling: the swap
    produced **0 new reds** over 1668 tests and 1252 constructions — it
    survives the extent guard (:math:`|\Gamma_+| = |\Gamma_-|` on every shipped
    quadrature) and changes no arithmetic. Asserted with ``is``, because
    :meth:`FunctionSpace.__eq__` compares ``(name, shape)`` and would accept a
    look-alike from the wrong face.
    """

    @pytest.mark.parametrize("case", _PERIODIC_CASES)
    def test_the_law_maps_the_PARTNERS_outflow_to_this_faces_inflow(
        self, case: str
    ) -> None:
        """``domain is Γ₊(partner)`` and ``codomain is Γ₋(face)``."""
        trace, realized, partner = _realized_periodic(case)
        assert partner == "xmax"
        assert realized.domain is trace.outflow_space(partner), (
            "periodic's domain is the PARTNER's outflow half — that is the "
            "whole content of the torus quotient, and the only off-diagonal "
            "block of B"
        )
        assert realized.codomain is trace.inflow_space("xmin")

    def test_the_binding_is_to_the_partner_not_to_THIS_face(self) -> None:
        r"""`[M]` :math:`\Gamma_+(\text{xmin})` and :math:`\Gamma_+(\text{xmax})`
        have IDENTICAL shape, so no shape or extent guard can tell them apart.

        This row is the load-bearing justification for binding spaces at all
        rather than checking lengths — and it is the catcher for a partner
        threading that quietly reverted to the installation face.
        """
        trace, realized, _ = _realized_periodic("gauss_legendre(8)")
        assert realized.domain is not trace.outflow_space("xmin")
        assert (
            trace.outflow_space("xmin").shape
            == trace.outflow_space("xmax").shape
        ), (
            "if these shapes ever differ, the cross-face claim stops being "
            "the space check's exclusive territory — re-read this gate"
        )

    def test_the_two_ends_are_DIFFERENT_spaces(self) -> None:
        """A wrap is an isomorphism between two spaces, never an endomorphism.

        The structural statement behind the ``IdentityOperator`` retirement,
        and it also rules out the both-ends-to-Γ₊ mis-binding, which the row
        above cannot see.
        """
        _, realized, _ = _realized_periodic("gauss_legendre(8)")
        assert realized.domain is not realized.codomain
        assert realized.domain is not None and realized.codomain is not None


class TestPeriodicChangedNoValues:
    r"""The migration-inheritance keystone for the wrap.

    The wrap's local permutation is ``arange`` (its linear part is the
    identity and the two index sets coincide — `[M]` an exact set equality on
    every shipped quadrature), so the uplift MUST be value-identical to the
    ``IdentityOperator`` it replaces or it changed physics rather than typing.
    """

    @pytest.mark.parametrize("case", _PERIODIC_CASES)
    def test_apply_and_transpose_are_the_identity_on_the_values(
        self, case: str
    ) -> None:
        """``array_equal`` — a gather is reduction depth 0, so bit-exact."""
        trace, realized, partner = _realized_periodic(case)
        rng = np.random.default_rng(20260807)
        x = rng.standard_normal(trace.outflow_space(partner).shape)
        np.testing.assert_array_equal(realized.apply(x), x)
        # The runtime narrowing IS the static permission (the operator
        # module's Design C). Asserting it here is also the claim that a
        # realized wrap HAS an honest transpose — PairedDeck declares it a
        # theorem, and its predecessor SpatialWrap answered False.
        assert adjointable(realized), (
            "the realized periodic law must expose a transpose; PairedDeck "
            "declares is_adjointable a THEOREM for a deck transformation"
        )
        y = rng.standard_normal(trace.inflow_space("xmin").shape)
        np.testing.assert_array_equal(realized.apply_transpose(y), y)

    def test_the_output_is_a_COPY_and_no_longer_ALIASES_the_input(
        self,
    ) -> None:
        r"""⚠ The one measured behavioural difference, recorded as a contract.

        `[M]` before step 7 the realized periodic law was
        ``IdentityOperator() & IdentityOperator()`` and ``apply(x) is x`` was
        **True**; a permutation returns a copy. Any consumer relying on
        identity-aliasing of the periodic face changes behaviour, and the
        campaign plan flagged it as a check-before-landing. Asserting it here
        turns "we checked once" into a standing contract.
        """
        trace, realized, partner = _realized_periodic("gauss_legendre(8)")
        x = np.arange(
            int(np.prod(trace.outflow_space(partner).shape)), dtype=float
        ).reshape(trace.outflow_space(partner).shape)
        out = realized.apply(x)
        assert out is not x, (
            "the realized wrap must return a COPY; an aliasing output means "
            "the IdentityOperator survived somewhere in the chain"
        )
        np.testing.assert_array_equal(out, x)


# ============================================================================
# G — ONE body
# ============================================================================


def test_both_deck_arms_enter_the_SAME_body(monkeypatch) -> None:
    r"""⭐ The Pattern-2 claim of the whole step, and nothing else can see it.

    A Mode-11 wrap counter: every *value* gate in this module stays green if a
    second private ``_periodic_kernel`` is re-added beside ``_deck_kernel``,
    because both would compute the same numbers. Only a counter on the shared
    body distinguishes "one source" from "two agreeing sources" — and the
    count is asserted **EXACTLY**, never ``> 0``, so a periodic arm that
    routed around the body would read 1 instead of 2.
    """
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    trace = _trace(quad, _X_FACES)
    realizer = SNBoundaryRealizer()

    calls: list[str] = []
    original = realizer_module._deck_kernel

    def counting(*args, **kwargs):
        calls.append(str(kwargs.get("law_key", "?")))
        return original(*args, **kwargs)

    monkeypatch.setattr(realizer_module, "_deck_kernel", counting)

    realizer.realize(
        ReflectiveBoundary(axis="x"),
        SNMethodSpace.for_face(quadrature=quad, face="xmin", trace=trace),
    )
    assert len(calls) == 1, (
        f"the specular arm must enter _deck_kernel exactly once; got {calls}"
    )
    realizer.realize(
        PeriodicBoundary(axis="x"),
        SNMethodSpace.for_face(quadrature=quad, face="xmax", trace=trace),
    )
    assert len(calls) == 2, (
        f"the PERIODIC arm did not enter _deck_kernel — it is still realized "
        f"by a body of its own, which is the twin path step 7 exists to "
        f"delete. Calls seen: {calls}"
    )


# ============================================================================
# The taxonomy reaches the realizer
# ============================================================================


def test_the_two_deck_types_realize_through_one_kernel_with_different_faces(
) -> None:
    r"""The self-paired / genuinely-paired split, observed at the realization.

    ``SelfPairedDeck.mirror`` consumes THIS face's :math:`\Gamma_+`;
    ``PairedDeck.wrap`` consumes the PARTNER's. One body, two faces — which is
    the step's outcome stated at the tier a consumer sees.
    """
    assert SelfPairedDeck.mirror(axis="x").domain_face("xmin") == "xmin"
    assert PairedDeck.wrap(axis="x").domain_face("xmin") == "xmax"

    trace, realized, partner = _realized_periodic("gauss_legendre(8)")
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    mirror_space = SNMethodSpace.for_face(
        quadrature=quad, face="xmin", trace=trace
    )
    mirrored = SNBoundaryRealizer().realize(
        ReflectiveBoundary(axis="x"), mirror_space
    )
    assert mirrored.domain is trace.outflow_space("xmin")
    assert realized.domain is trace.outflow_space(partner)
    assert mirrored.codomain is realized.codomain, (
        "both laws deliver into the SAME half-trace — Γ₋(xmin) — and differ "
        "only in which Γ₊ they consume. That is the whole taxonomy."
    )
