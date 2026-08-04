r"""Bit-exactness of the realized :math:`O(3)` operators (campaign step G4).

A distinct claim class from :mod:`tests.numerics.test_symmetry`, which
gates the *group structure* — which subgroup contains which, and which
rule is invariant under which group. This file gates the **floating-point
realization** of the same operators: how exactly the matrices are built,
and how exactly they carry a real quadrature's node set onto itself.

Why the two are separate questions, and why this one needs its own gates:
orthogonality, determinant, closure and group order are all preserved by a
*less accurate* construction of the same group, so nothing in the
structural file can see an operator set drifting into a wider tolerance.
The drift is silent, and it is absorbed by the very node-match window the
checker uses to decide invariance.

The construction under test
---------------------------

Both azimuthal families are built from
:func:`~orpheus.numerics.roots_of_unity.roots_of_unity` — the SAME
generator the azimuthal *rules* are built from — rather than from
:func:`numpy.cos` of an angle:

* :math:`C_n` from the :math:`n`-th roots as exact circle points;
* the vertical mirrors as the coset :math:`C_n\,\sigma_0`, so the mirror
  set is the rotation set carried by ONE bit-exact coordinate reflection.

What is achievable, measured
----------------------------

Exactness is **not** available at every order, and the reason is
structural rather than a defect to be fixed later. Carrying the node at
azimuth :math:`2\pi j/n` to the node at :math:`2\pi(j+k)/n` requires

.. math::

    \cos a \cos b - \sin a \sin b = \cos(a + b)

to hold in IEEE-754, and the angle-addition formula is not a
floating-point theorem. So the landing residual is exactly zero only when
every root involved is an axis point — :math:`n_\varphi \in \{1, 2, 4\}` —
and is a few ULP otherwise. `[M]` over :math:`n_\varphi = 1 \ldots 24`:
zero at 1/2/4, and between ``1.1e-16`` and ``3.3e-16`` elsewhere.

The campaign plan's acceptance line for this step read
":math:`D_{nh}(n_\varphi)` exact on BOTH sides". That is achievable
literally at the quarter-turn orders and provably *unachievable* at the
rest; the honest criterion is the pair below — exact where exactness
exists, and the mirror half no longer worse than the rotation half
anywhere.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.quadrature import product_mu_phi
from orpheus.numerics.symmetry import _cyclic_ops, _vertical_mirrors

pytestmark = pytest.mark.foundation


#: Values a symmetric-angle matrix entry is allowed to be *exactly*.
#: ``sqrt(1/2)`` is the 45-degree diagonal, shared by both components.
_SQRT_HALF = float(np.sqrt(0.5))
_EXACTLY_REPRESENTABLE = (0.0, 1.0, -1.0, _SQRT_HALF, -_SQRT_HALF)

#: The orders whose every root of unity is an axis point, hence the only
#: ones where the angle-addition step is exact. See the module docstring.
QUARTER_TURN_ORDERS = (1, 2, 4)


def _landing_residual(nodes: np.ndarray, ops) -> float:
    r"""``max_g max_i min_j |g(x_i) - x_j|`` — how far the operator set
    misses carrying the node set onto itself.

    Deliberately the raw distance rather than
    :meth:`~orpheus.geometry.transformation.RigidMotion.permutes`, whose
    ``atol`` would absorb exactly the quantity being measured.
    """
    worst = 0.0
    for g in ops:
        moved = g.on_points(nodes)
        distance = np.linalg.norm(moved[:, None, :] - nodes[None, :, :], axis=2)
        worst = max(worst, float(np.max(np.min(distance, axis=1))))
    return worst


def _sphere_nodes(n_mu: int, n_phi: int) -> np.ndarray:
    measure, _ = product_mu_phi(n_mu, n_phi)
    return np.asarray(measure.nodes, dtype=float)


# ============================================================================
# The exactness that IS achievable
# ============================================================================


@pytest.mark.parametrize("n_mu", [2, 4, 6])
@pytest.mark.parametrize("n_phi", QUARTER_TURN_ORDERS)
def test_quarter_turn_orders_carry_the_rule_EXACTLY(n_mu: int, n_phi: int) -> None:
    r"""At :math:`n_\varphi \in \{1,2,4\}` the residual is ``0.0``, not ``~0``.

    Every root of unity involved is an axis point, so the rotation matrix
    entries are exactly :math:`0` and :math:`\pm 1` and carrying a node is
    a coordinate permutation with sign flips — no arithmetic to round.
    Before G4 the same orders measured ``2.1e-16`` (rotations) and
    ``4.7e-16`` (mirrors), because the operators were built from
    :func:`numpy.cos` of a rounded :math:`\pi/2`.

    Asserted at three polar orders, since the azimuthal claim must not
    depend on the polar factor.
    """
    nodes = _sphere_nodes(n_mu, n_phi)
    rotations = _landing_residual(nodes, _cyclic_ops(n_phi))
    mirrors = _landing_residual(nodes, _vertical_mirrors(n_phi))
    assert rotations == 0.0, (
        f"C_{n_phi} misses product({n_mu},{n_phi}) by {rotations:.3e}; at a "
        f"quarter-turn order every operator entry is exactly 0 or +-1, so "
        f"any non-zero residual means the operators are no longer built "
        f"from exact circle points"
    )
    assert mirrors == 0.0, (
        f"the vertical mirrors miss product({n_mu},{n_phi}) by "
        f"{mirrors:.3e}; they are the coset C_n . sigma_0 with sigma_0 a "
        f"bit-exact coordinate reflection, so they cannot be less exact "
        f"than the rotations"
    )


def test_the_residual_instrument_can_actually_read_non_zero() -> None:
    r"""Positive control: SOME order in the sweep measures non-zero.

    Without this leg the exactness gates above are unfalsifiable by a
    broken measurement — a ``_landing_residual`` that returned ``0.0``
    unconditionally would make every one of them pass. (`vv-principles`
    anti-pattern #17: verify the instrument on a known-positive before
    trusting any negative it reports.)

    This pins the INSTRUMENT, not a limitation: it asserts that the sweep
    contains a non-zero, never that any particular order must stay
    inexact. Making an order exact is an improvement and reds nothing
    here as long as one non-zero remains.
    """
    residuals = {
        n_phi: _landing_residual(_sphere_nodes(4, n_phi), _cyclic_ops(n_phi))
        for n_phi in range(1, 13)
    }
    assert any(r > 0.0 for r in residuals.values()), (
        "every order measured exactly 0.0 — the residual instrument is "
        f"reading a constant, not a distance. Measured: {residuals}"
    )


# ============================================================================
# The coset consequence — the mirror half is no longer the worse half
# ============================================================================


@pytest.mark.parametrize("n_phi", list(range(1, 17)))
def test_mirrors_are_exactly_as_accurate_as_the_rotations(n_phi: int) -> None:
    r"""``residual(sigma_v set) == residual(C_n set)``, bit-for-bit.

    The consequence of building the mirrors as :math:`C_n\,\sigma_0`:
    :math:`\sigma_0` is a coordinate mirror, hence an exactly-representable
    signed diagonal, so composing with it is a column sign flip that
    introduces no round-off — and the rule's own node set is carried by
    :math:`\sigma_0` exactly (that is the mirror guarantee
    :func:`~orpheus.numerics.roots_of_unity.roots_of_unity` exists to
    provide). The mirror residual is therefore the rotation residual,
    reindexed.

    This is the falsifiable half of the change, and it is what the
    superseded construction could not give: building each mirror from a
    normal at the HALF angle :math:`k\pi/n + \pi/2` used a root of unity
    of order :math:`4n` — a different generator from the rule's own — and
    then normalised it, leaving the mirrors 2--5x LESS accurate than the
    rotations (`[M]` ``1.26e-15`` vs ``7.1e-16`` at
    :math:`n_\varphi = 16`).

    Equality, not a tolerance: two floats computed from the same
    quantities either agree exactly or the construction changed.
    """
    nodes = _sphere_nodes(4, n_phi)
    rotations = _landing_residual(nodes, _cyclic_ops(n_phi))
    mirrors = _landing_residual(nodes, _vertical_mirrors(n_phi))
    assert mirrors == rotations, (
        f"n_phi={n_phi}: mirrors miss by {mirrors:.6e} but rotations by "
        f"{rotations:.6e}. The vertical mirrors are supposed to BE the "
        f"rotations composed with one exact coordinate reflection; a "
        f"difference means either the coset construction was replaced by "
        f"an independent one, or the rule is no longer exactly closed "
        f"under the phi=0 mirror"
    )


# ============================================================================
# Matrix-level exactness
# ============================================================================


@pytest.mark.parametrize("n_phi", QUARTER_TURN_ORDERS)
def test_quarter_turn_operators_are_exact_signed_permutations(n_phi: int) -> None:
    r"""Every entry of every :math:`C_n` and :math:`\sigma_v` matrix at
    :math:`n_\varphi \in \{1,2,4\}` is exactly :math:`0` or :math:`\pm 1`.

    The matrix-level statement behind the landing gate, and the one that
    localises a regression: a residual can drift for either the operator
    or the rule, but this reads the operator alone. `[M]` before G4,
    :math:`C_4` had 30 of 36 entries exact and the mirrors 26 of 36 --
    ``np.cos(np.pi/2)`` is ``6.1e-17``, not ``0``.
    """
    for family, ops in (("C_n", _cyclic_ops(n_phi)),
                        ("sigma_v", _vertical_mirrors(n_phi))):
        entries = np.concatenate([np.asarray(g.linear).ravel() for g in ops])
        inexact = entries[~np.isin(entries, (0.0, 1.0, -1.0))]
        assert inexact.size == 0, (
            f"{family} at n_phi={n_phi} carries {inexact.size} entries that "
            f"are not exactly 0 or +-1: {np.unique(inexact)}"
        )


def test_eighth_turn_operators_hit_the_diagonal_constant_exactly() -> None:
    r"""At :math:`n_\varphi = 8` every entry is exactly :math:`0`,
    :math:`\pm 1`, or :math:`\pm\sqrt{1/2}`.

    The 45-degree diagonal is the octant fold's fixed point, where
    :math:`\cos` and :math:`\sin` must be the SAME float — which
    ``np.cos(np.pi/4)`` and ``np.sin(np.pi/4)`` are not (they differ by
    1 ULP). `[M]` 68 of 72 entries exact per family after G4, against 54
    (rotations) and 47 (mirrors) before; the 4 residual entries are the
    genuine eighth-turn cosines, which have no exact representation.
    """
    for family, ops in (("C_8", _cyclic_ops(8)),
                        ("sigma_v(8)", _vertical_mirrors(8))):
        entries = np.concatenate([np.asarray(g.linear).ravel() for g in ops])
        exact = int(np.isin(entries, _EXACTLY_REPRESENTABLE).sum())
        assert exact == 68, (
            f"{family}: {exact}/72 entries exactly representable, expected "
            f"68. Fewer means the operators left the roots-of-unity "
            f"construction; more means an exactness improvement that "
            f"should be recorded here rather than silently absorbed"
        )
