r"""The deck arrow is the PREIMAGE of the deck's motion: the gather convention of ``ScaledOrdinatePermutation`` (step G6.3-7, re-posed at wave 2 S2b).

Sibling of ``test_specular_deck_chain.py`` (the mirror's binding) and ``test_lambertian_chain.py`` (the two-link response chain).
This file gated ``_deck_kernel``, the one body every deck transformation's trace arrow realized through, which derived the arrow
from the deck element's **rigid motion**. The carve retired that builder (decision 45: the typed operators OWN the one body of each
builder): the arrow is now
:class:`~orpheus.sn.boundary.realized_operators.ScaledOrdinatePermutation`, whose constructor takes the ordinate permutation
``pi`` (the caller derives it: ``ordinate_pairing`` for a mirror, the identity for a translation wrap) and whose operator gathers
``pi^{-1}[inflow]`` at the ordinate's position inside :math:`\Gamma_+(f')` (``_crossing``).

``h`` is the deck element carrying the DOMAIN face onto the INSTALLATION face, :math:`h(f') = f`. Invariance of :math:`\psi` under
the deck group gives

.. math::

   \gamma_-\psi|_f(\Omega) \;=\; \gamma_+\psi|_{f'}(h^{-1}\Omega),

so with :meth:`~orpheus.geometry.transformation.RigidMotion.permutes`' convention :math:`h(x_i) = x_{\pi(i)}` the row for inflow
ordinate ``j`` gathers global ordinate :math:`\pi^{-1}(j)`, remapped into :math:`\Gamma_+(f')`'s local index by ``to_local``.

**What this file keeps, and why only this.** Every SHIPPED deck is an involution (the mirror) or the identity (the wrap), for which
:math:`\pi = \pi^{-1}`: the two gather conventions coincide there, and only a non-involutive pairing (a sector rotation, the
latent #178 consumer) can tell a gather of :math:`\pi^{-1}` from one of :math:`\pi`. The rows below pose the convention on the six
admissible 90-degree rotation fixtures, where the wrong convention is *unrealizable* rather than silently wrong, so the latent
consumer inherits it correct. ⭐ Three measured facts shaped them (`[M]` 2026-08-07, kept):

**1.** The wrong gather direction is REFUSED, not silently wrong: on all six rotation fixtures :math:`\pi^{-1}[\text{inflow}]` lands
inside :math:`\Gamma_+(f')` and :math:`\pi[\text{inflow}]` lands outside it, so ``to_local`` raises (structural: :math:`h(\Gamma_-(f)) =
h^2(\Gamma_+(f'))`, which equals :math:`\Gamma_+(f')` iff h is an involution). Consequence for credit: a mutation that reds by RAISING
attributes nothing to the value rows; the discrimination is pinned by :meth:`TestTheGatherIsThePreimage.test_the_forward_gather_is_OUT_of_range`.

**2.** ``product(4,4)`` is the DEGENERATE rotation fixture: its local permutation is exactly ``arange``, like the wrap's, so a
keystone written there is green under ``local := np.arange(n)``. The keystone runs on ``product(4,8)``; ``product(4,4)`` is kept as a
labelled degeneracy control.

**3.** Tolerances per law: the witness residual :math:`\max|x_g - h^{-1}x_j|` is **exactly 0.0** for the mirror (a signed permutation is
exact in IEEE) and the wrap, and **<= 1.11e-16** for the 90-degree rotation (``cos(pi/2)`` round-off), so the mirror and wrap rows
assert ``array_equal`` and only the rotation rows carry ``atol = 1e-15``.

**Retired with the builder** (demotion check): ``TestThePermutationIsTheLinearParts`` (the builder took a motion and had to read its
LINEAR part; the realizer now never hands a translation motion to ``ordinate_permutation``: a wrap's permutation is the identity by
construction); ``TestTheDeckKernelRefuses`` and the duplicated-node row (refusals of the motion-to-permutation matcher, owned by
``Quadrature.ordinate_permutation`` and ``ordinate_pairing``, gated in ``tests/gates/numerics/test_ordinate_pairing.py`` and the
quadrature tests; the table's refusal of an uncertified mirror is ``test_sn_boundary_realizer.py``); ``TestPeriodicIsBoundToThePartner`` /
``TestPeriodicChangedNoValues`` (re-posed in ``test_sn_boundary_realizer.py::TestRealizeTranslation`` and the METRIC class);
``test_both_deck_arms_enter_the_SAME_body`` (one class now: unspellable).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry.transformation import RigidMotion
from orpheus.numerics.face_layout import FaceLayout
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.spaces import AngularTraceSpace
from orpheus.sn.boundary.realized_operators import ScaledOrdinatePermutation
from orpheus.sn.mesh.method_space import SNMethodSpace

pytestmark = pytest.mark.l1

_NG = 2

#: The measured witness tolerance for a rotation deck (see the module docstring). The mirror and the wrap are asserted bit-exact.
_ROTATION_ATOL = 1e-15

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
        [(f, (int(quad.n_ordinates), _NG)) for f in faces]
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
    gamma_out = trace.outflow_space(domain_face)
    return quad, trace, method_space, gamma_out, quad.measure.support.ambient_points(quad.measure.nodes)


def _kernel(quadrature, method_space, gamma_out, *, motion, domain_face, law_key):
    """The PermutationOperator of the deck arrow ``Gamma_+(domain_face) -> Gamma_-(method_space.face)`` of ``motion``.

    THE ONE CALL SITE: the permutation is the one the motion induces on the ORDINATES (the linear part's action:
    a translation moves no direction), and ``ScaledOrdinatePermutation`` gathers its inverse. ``gamma_out`` and ``law_key`` are
    carried only so the rows read as they always did.
    """
    permutation = quadrature.ordinate_permutation(RigidMotion(motion.linear))
    op = ScaledOrdinatePermutation(alpha=1.0, permutation=permutation, face=method_space.face, partner=domain_face, trace=method_space.trace)
    return op.inner.ops[0]


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
    return np.asarray(gamma_out.ordinate_indices, dtype=np.intp)[np.asarray(kernel.perm)]


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
            motion=RigidMotion.translation_by([1.0, 0.0, 0.0]),
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
        nodes = quad.measure.support.ambient_points(quad.measure.nodes)
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
            np.sort(np.asarray(gamma_out.ordinate_indices)), gather
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
            gamma_out.to_local(forward)
        # …and the correct one is accepted, so the refusal above is about the
        # DIRECTION and not about the space refusing everything.
        gamma_out.to_local(backward)


def _parse_product(case: str) -> tuple[int, int]:
    """``"product(4,8)" -> (4, 8)``."""
    inner = case[case.index("(") + 1: case.index(")")]
    n_mu, n_phi = (int(t) for t in inner.split(","))
    return n_mu, n_phi

