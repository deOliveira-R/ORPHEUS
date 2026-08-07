r"""The specular RESPONSE, pinned against the geometric deck transformation.

G6.3 step 4 (#330), from a user ruling of 2026-08-04. Two statements, and the
second is a gate the tree did not have.

**A deck transformation admits NO attenuation — and that is the PROOF it is
geometry.** :class:`~orpheus.geometry.boundary.ReflectiveBoundary` is a symmetry
plane: a quotient of the domain, a theorem imposed, carrying zero physics. A
theorem cannot absorb. So an attenuated "reflective" law is not a dimmer mirror
— it is an :class:`~orpheus.geometry.boundary.AlbedoBoundary` with a
:class:`~orpheus.geometry.boundary.SpecularReturn` closure, wearing the geometry
costume. ``_factors.py`` already flags that row as a deliberate violation; this
module supplies the *consequence* form of the argument: **if it can be
attenuated, it was never a quotient.**

⭐ **And the unattenuated specular RESPONSE must equal the geometric deck
transformation.** The two are different KINDS of object — one is a constitutive
stand-in for a real mirror, the other a theorem about the domain — and they
coincide at exactly one point, :math:`\alpha = 1`. That coincidence is
checkable, and checking it is what licenses the response as a stand-in at all.

**Why this is a real cross-check and not a tautology.** The two DERIVATIONS
are genuinely independent:

* the **geometric** side applies
  :meth:`~orpheus.geometry.boundary.SelfPairedDeck.mirror`'s
  :class:`~orpheus.geometry.transformation.RigidMotion` — the ``G1``–``G5``
  core, verified against pure math — to the quadrature's direction cosines and
  reads off the induced permutation with this module's OWN
  argmin + exactness + injectivity discipline (no call into ``preserves``);
* the **response** side is ``quadrature.ordinate_permutation(motion)``, the
  quadrature's own certified pairing derivation — the one source production
  realization also reads (G6.3 step 7; until §7d.3 this side was the
  construction-time ``reflection_index`` table).

⚠ **Scope, narrowed at §7d.3 — independence has TWO axes, and this file
holds one.** Both sides consume ONE ``SelfPairedDeck.mirror(axis).motion``
object, so the INPUT — the axis-letter → mirror-normal resolution — is
shared, and this file no longer cross-checks that convention: `[M]`
swapping x↔y inside ``_mirror_motion`` leaves all rows here GREEN while
78 sibling gates red. The letter-convention cross-check lives in the
reference helper's deliberately-local axis map
(``tests/_harness/references.py`` ``_AXIS_INDEX``) and the gates it
feeds. What THIS file pins is the DERIVATION: ``preserves``'s match
machinery against an independently-written match, on the same motion.

`[M]` the two derivations agree EXACTLY (``atol=1e-12``, and the
permutation is a genuine bijection) on ``gauss_legendre(4)``,
``gauss_legendre(8)``, ``product(4,4)`` (x and y), ``lebedev(17)`` and
``level_symmetric(6)`` (z).

Beyond :math:`\alpha = 1` the two part company for good, and the user's framing
is worth keeping: **even a scalar attenuation is already an imposed response** —
a real mirror's attenuation is rarely uniform in angle. The scalar is a
modelling choice, not a measurement, and the moment it is applied there is no
geometry left to check against.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry.boundary import ReflectiveBoundary
from orpheus.geometry.boundary._factors import SelfPairedDeck
from orpheus.numerics.quadrature import Quadrature

pytestmark = pytest.mark.l1

_CASES = {
    "gauss_legendre(4) x": (lambda: Quadrature.gauss_legendre(n_ordinates=4), "x"),
    "gauss_legendre(8) x": (lambda: Quadrature.gauss_legendre(n_ordinates=8), "x"),
    "product(4,4) x": (lambda: Quadrature.product(n_mu=4, n_phi=4), "x"),
    "product(4,4) y": (lambda: Quadrature.product(n_mu=4, n_phi=4), "y"),
    "lebedev(17) x": (lambda: Quadrature.lebedev(order=17), "x"),
    "lebedev(17) y": (lambda: Quadrature.lebedev(order=17), "y"),
    "level_symmetric(6) z": (lambda: Quadrature.level_symmetric(sn_order=6), "z"),
}


def _directions(quad) -> np.ndarray:
    """``(N, 3)`` direction cosines — uniform across 1-D and 3-D rules.

    A 1-D slab rule stores bare :math:`\\mu` and has no ``mu_y`` / ``mu_z``; the
    missing components are genuinely zero, so padding is the honest lift rather
    than a convenience. (Reading ``quad.nodes`` instead does NOT work: it is
    ``(N,)`` for a 1-D rule and ``(3, N)`` for a 3-D one, and treating the
    former as a coordinate array silently makes a 4-dimensional ambient space.)
    """
    return np.column_stack([
        np.zeros(quad.N) if getattr(quad, a, None) is None
        else np.asarray(getattr(quad, a), dtype=float)
        for a in ("mu_x", "mu_y", "mu_z")
    ])


def _geometric_permutation(quad, axis: str) -> tuple[np.ndarray, np.ndarray]:
    """The permutation the DECK TRANSFORMATION induces on the ordinates.

    Applies the rigid motion and matches each image back to a node. Returns
    ``(permutation, residuals)`` — the residuals so the caller can assert the
    match is EXACT rather than merely nearest, which is the difference between
    a permutation and a relation (``vv-principles`` anti-pattern #14).
    """
    nodes = _directions(quad)
    motion = SelfPairedDeck.mirror(axis=axis, dimension=3).motion
    moved = nodes @ np.asarray(motion.linear).T
    distances = np.linalg.norm(moved[:, None, :] - nodes[None, :, :], axis=-1)
    induced = np.argmin(distances, axis=1)
    return induced, distances[np.arange(len(nodes)), induced]


@pytest.mark.parametrize("case", list(_CASES))
def test_the_deck_transformation_induces_the_quadratures_own_table(case):
    r"""⭐ ``σ_axis`` applied to the nodes == ``quadrature.ordinate_permutation``.

    The pin: the constitutive stand-in and the geometric theorem agree at
    :math:`\alpha = 1`, by two derivations that never consult each other.

    Asserted with the match's EXACTNESS and BIJECTIVITY, not just the equality:
    a nearest-neighbour loop computes a *relation*, and a many-to-one map would
    satisfy the equality check while being no permutation at all (ERR-073's
    shape). Both legs are cheap and both are necessary.

    ⭐ **The exactness leg is not belt-and-braces — it catches a class the
    equality leg cannot see.** `[M]` a SCALED mirror ``diag(-1,1,1)·1.1``
    induces the CORRECT permutation (nearest neighbour is still the mirror
    partner) while its images are not nodes at all: ``matches_table=True``,
    ``exact=False``. Equality alone would pass it. What the residual catches is
    exactly the **measure-preserving** half of the deck-transformation
    definition — a scaling is not an isometry, so it is not a deck element
    however right its combinatorics look. Mutations checked: wrong axis,
    inversion, identity, and that scaling; the first three fail on equality,
    the fourth on exactness alone.
    """
    quad = _CASES[case][0]()
    axis = _CASES[case][1]
    induced, residuals = _geometric_permutation(quad, axis)

    assert np.allclose(residuals, 0.0, atol=1e-12), (
        f"the mirror image of a node is not a node (max residual "
        f"{residuals.max():.3e}) — the quadrature is not symmetric under "
        f"σ_{axis}, so no permutation exists to compare"
    )
    assert len(set(induced.tolist())) == len(induced), (
        "the induced map is not injective — it is a relation, not a permutation"
    )
    motion = SelfPairedDeck.mirror(axis=axis, dimension=3).motion
    pi = quad.ordinate_permutation(motion)
    assert pi is not None, (
        f"production's pairing derivation refuses σ_{axis} on a rule the "
        f"geometric side just proved symmetric under it"
    )
    np.testing.assert_array_equal(induced, pi.indices)


@pytest.mark.parametrize("case", list(_CASES))
def test_the_geometric_permutation_is_an_involution(case):
    r"""``σ² = id`` on the ordinates — the theorem, read on the discretisation.

    A mirror is an involution in the continuum; the discrete table inherits it
    only if the rule is symmetric. This is the same claim ERR-044 guards at the
    law tier, asserted here on the GEOMETRIC side so the two tiers are pinned
    independently rather than through one shared table.
    """
    quad = _CASES[case][0]()
    induced, _ = _geometric_permutation(quad, _CASES[case][1])
    np.testing.assert_array_equal(induced[induced], np.arange(quad.N))


def test_an_attenuated_reflective_law_is_NOT_a_deck_transformation():
    r"""⭐ A deck transformation admits no attenuation — the consequence form.

    ``ReflectiveBoundary`` still accepts an ``albedo``, which ``_factors.py``
    flags as a deliberate, visible violation of "exactly one of ``G``, ``R`` is
    non-trivial": at :math:`\alpha \neq 1` BOTH factors are non-trivial and the
    object is an ``AlbedoBoundary(α, SpecularReturn(axis))`` in a geometry
    costume.

    This gate states the same thing as a **consequence** rather than a
    taxonomy: a quotient of the domain is a theorem, and a theorem cannot
    absorb — so attenuability is itself the discriminator. It is pinned rather
    than left as prose because the parameter is still reachable in Python (it is
    unreachable from a ``BC(...)`` tag, which hard-codes :math:`\alpha = 1`),
    and its retirement is deferred to phase B5. When that lands, this gate
    should red and be replaced by the constructor's refusal.
    """
    unattenuated = ReflectiveBoundary(axis="x", albedo=1.0)
    attenuated = ReflectiveBoundary(axis="x", albedo=0.7)

    # The geometry factor is the SAME mirror in both — attenuation is not a
    # property of the deck element, which is exactly the problem.
    assert unattenuated.geometry_map == attenuated.geometry_map

    # Only the unattenuated one is a pure quotient: at α=1 the response is
    # trivial, so the law's entire content is the deck transformation.
    assert unattenuated.response_kernel.amplitude == 1.0
    assert attenuated.response_kernel.amplitude == 0.7, (
        "the attenuation surfaced somewhere other than the response kernel — "
        "if it is not in R it has been smuggled into the geometry, which is "
        "the conflation this gate exists to refuse"
    )
