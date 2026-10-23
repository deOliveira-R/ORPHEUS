r"""Equivariance certification of a quadrature under a rigid motion.

:func:`ordinate_pairing` answers one question, keyed on
``(quadrature, motion)`` and nothing else (ruling 7 of ``layer_dag.md``): *is
the ordinate permutation this motion induces a certified pairing?* It carries
no law argument and no boundary vocabulary, so it lives in the mathematics
layer and every consumer that needs a certified mirror pairing (a mirror deck,
a specular chain) reads this one array.

Design
======

The function FOLDS onto :meth:`~orpheus.numerics.quadrature.Quadrature.ordinate_permutation`,
the single source of "which permutation does a motion induce":

* ``None`` (no bijective, weight-preserving match, ERR-073) becomes
  :class:`QuadratureNotClosedUnderMotionError`;
* on the RETURNED permutation three claims are asserted independently, because
  the matcher's built-in certificate (bijection plus bare-weight equality)
  does not subsume them:

  - **measure** (ERR-042): the trace measure :math:`m_n = w_n\,|\mu_{a,n}|`
    satisfies :math:`m_{\pi(n)} = m_n` to relative tolerance
    :data:`_MEASURE_RTOL`. The matcher compares the bare :math:`w`, not the
    cosine-weighted one;
  - **involution** (ERR-044): :math:`\pi\circ\pi = \mathrm{id}`. A bijection
    need not be its own inverse;
  - **crossing** (ERR-045): every non-tangential ordinate pairs with the
    opposite sign class on the axis :math:`a`.

The three are independent: each passes on a table the other two reject (the
catalogue's lesson, ERR-045), which is why none is dropped for another.

The mirror plane is stated ONCE, by the motion: the axis :math:`a` is the one
whose diagonal entry of ``motion.linear`` is :math:`-1` (all other entries
:math:`+1`). A motion that is not an axis reflection (rotation, translation,
identity, oblique reflection, inversion) names no mirror plane and raises the
base :class:`OrdinatePairingError`. There is no ``axis`` argument, so a
disagreeing pair (a reflection about x with ``axis="y"``) is unspellable.

The trace measure reads
:meth:`~orpheus.numerics.quadrature.Quadrature.mean_axis_cosine`, not
``axis_cosines``: this is a flux question, so on an axis a 1-D rule has
suppressed, zero is the answer.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from orpheus.geometry.faces import AXIS_NAMES
from orpheus.numerics.tangential import TANGENTIAL_EPS

if TYPE_CHECKING:
    from orpheus.geometry.transformation import Permutation, RigidMotion

    from .directional import Quadrature

__all__ = [
    "OrdinatePairingError",
    "PairingDoesNotCrossError",
    "PairingNotInvolutiveError",
    "PairingNotMeasurePreservingError",
    "QuadratureNotClosedUnderMotionError",
    "ordinate_pairing",
]

#: Relative tolerance of the ERR-042 measure check. Partners come from the
#: same symmetric node/weight sets, so a correct table agrees to a few ULP and
#: a mispaired weight class violates at O(1).
_MEASURE_RTOL: float = 1e-12


class OrdinatePairingError(ValueError):
    """The motion has no certified ordinate pairing on the quadrature."""


class QuadratureNotClosedUnderMotionError(OrdinatePairingError):
    """No bijective weight-preserving match of the ordinates onto their
    images (e.g. an odd-``n_phi`` product rule has no x-mirror closure)."""


class PairingNotMeasurePreservingError(OrdinatePairingError):
    """ERR-042: the pairing does not preserve :math:`w\\,|\\mu_a|`."""


class PairingNotInvolutiveError(OrdinatePairingError):
    """ERR-044: the pairing is not its own inverse."""


class PairingDoesNotCrossError(OrdinatePairingError):
    """ERR-045: a non-tangential ordinate is not paired with the opposite
    sign class on the motion's axis."""


def _mirror_axis(motion: "RigidMotion") -> int:
    r"""The axis index of the mirror plane the motion states, or a refusal.

    The linear part must be :math:`\mathrm{diag}(1,\dots,-1,\dots,1)` with
    exactly one :math:`-1`.
    """
    m3 = motion if motion.dimension == 3 else motion.embedded_in(3)
    diagonal = np.diag(m3.linear)
    if not np.allclose(m3.linear, np.diag(diagonal), atol=1e-12, rtol=0.0):
        raise OrdinatePairingError(
            "the motion is not an axis reflection: its linear part is not "
            "diagonal, so it names no mirror plane normal to a coordinate axis"
        )
    negative = np.isclose(diagonal, -1.0, atol=1e-12, rtol=0.0)
    positive = np.isclose(diagonal, 1.0, atol=1e-12, rtol=0.0)
    if int(negative.sum()) != 1 or not (negative | positive).all():
        raise OrdinatePairingError(
            f"the motion is not an axis reflection: its linear diagonal "
            f"{diagonal.tolist()} does not hold exactly one -1 and the rest +1"
        )
    return int(np.flatnonzero(negative)[0])


def ordinate_pairing(quadrature: "Quadrature", motion: "RigidMotion") -> "Permutation":
    r"""The certified ordinate permutation of an axis reflection.

    Parameters
    ----------
    quadrature:
        The rule. Only ``n_ordinates``, ``weights``, ``measure``,
        ``axis_cosines``, ``mean_axis_cosine`` and ``ordinate_permutation``
        are read.
    motion:
        An axis reflection; the mirror axis is the one it names.

    Returns
    -------
    Permutation
        :math:`\pi` with ``pi.indices[i]`` the index of ordinate :math:`i`'s
        image, the very array
        :meth:`~orpheus.numerics.quadrature.Quadrature.ordinate_permutation`
        returns, after the three checks.

    Raises
    ------
    OrdinatePairingError
        The motion is not an axis reflection.
    QuadratureNotClosedUnderMotionError
        The rule is not closed under the motion.
    PairingNotMeasurePreservingError
        ERR-042.
    PairingNotInvolutiveError
        ERR-044.
    PairingDoesNotCrossError
        ERR-045.
    """
    axis_index = _mirror_axis(motion)
    axis = AXIS_NAMES[axis_index]
    permutation = quadrature.ordinate_permutation(motion)
    if permutation is None:
        raise QuadratureNotClosedUnderMotionError(
            f"no ordinate pairing about axis {axis!r} exists on this "
            f"quadrature: the mirror does not permute its weighted ordinate "
            f"set (no bijective, weight-preserving match, e.g. an odd-n_phi "
            f"product rule has no x-mirror closure)"
        )
    pi = permutation.indices
    mu = quadrature.mean_axis_cosine(axis_index)

    trace_measure = quadrature.weights * np.abs(mu)
    partner_measure = trace_measure[pi]
    if not np.allclose(partner_measure, trace_measure, rtol=_MEASURE_RTOL, atol=0.0):
        worst = int(np.argmax(np.abs(partner_measure - trace_measure)))
        raise PairingNotMeasurePreservingError(
            f"the pairing about axis {axis!r} does not preserve the "
            f"direction-cosine measure w*|mu_{axis}|: ordinate {worst} carries "
            f"m={trace_measure[worst]:.6e} but its partner {int(pi[worst])} "
            f"carries m={partner_measure[worst]:.6e} (ERR-042)"
        )

    if not np.array_equal(pi[pi], np.arange(int(quadrature.n_ordinates))):
        raise PairingNotInvolutiveError(
            f"the pairing about axis {axis!r} is not an involution (ERR-044)"
        )

    partner_mu = mu[pi]
    active = np.abs(mu) > TANGENTIAL_EPS
    partner_opposite = (np.abs(partner_mu) > TANGENTIAL_EPS) & (
        np.sign(partner_mu) != np.sign(mu)
    )
    bad = active & ~partner_opposite
    if bad.any():
        worst = int(np.flatnonzero(bad)[0])
        raise PairingDoesNotCrossError(
            f"the pairing about axis {axis!r} maps ordinate {worst} "
            f"(mu_{axis}={mu[worst]:+.6e}) to ordinate {int(pi[worst])} "
            f"(mu_{axis}={partner_mu[worst]:+.6e}), the same sign class "
            f"instead of the outflow partner (ERR-045)"
        )
    return permutation
