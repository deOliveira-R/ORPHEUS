r"""The gradings of the characteristic reference: where a quadrature puts its piece ends, and how many nodes a piece takes.

Four laws, each the one place it is spelled:

* **geometric** (:func:`graded_ends`): layers at :math:`w\rho^j` toward an
  end, the panels toward walls and interfaces and the halvings of the line
  rules;
* **exponential** (:func:`exponential_ends`): ends at a first depth times
  :math:`\rho^k` mean free paths toward where an attenuation concentrates,
  until it vanishes (:data:`VANISHING_DEPTH`), over a family of attenuations
  at once, at the :data:`DYADIC` first depth and ratio. It grades the rules
  over lines (the polar angle and the cosine toward the normal, toward
  grazing; :mod:`.lines`), coordinates in which the attenuation's exponent
  is not linear. Along a line the attenuation cuts nothing: it is the
  measure of each piece's rule (:mod:`.attenuation_rule`);
* **hp**: toward a singularity off the interval, until each piece is no
  wider than its distance to it over a tightness, so a Gauss rule converges
  geometrically on every piece. :func:`halvings` halves toward the impact
  rule's next radius and :math:`b = 0`, and the direction rules toward
  grazing; :func:`branch_walk` walks each slot of a line toward the orbit
  coordinate's complex branch points, at the tightness :data:`TIGHTNESS`;
* **the node count** (:func:`gauss_nodes`): the least number of nodes of a
  piece's Gauss rule that holds :data:`ALONG_LINE_ERROR`, from the Bernstein
  ellipse through the piece's nearest singularity.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

from .attenuation_rule import BAND

#: The optical depth past which an attenuation vanishes in double precision (e^-64): the last graded end reaches it.
VANISHING_DEPTH = 64.0


@dataclass(frozen=True)
class ExponentialGrading:
    r"""Where an exponential grading puts its ends, in mean free paths of the attenuation: a first depth and a ratio (:func:`exponential_ends`)."""

    first: float
    ratio: float


#: The dyadic grading, first depth 1 and ratio 2: for a coordinate in which the attenuation's exponent is not linear (the
#: polar angle and the cosine toward the normal, :mod:`.lines`), where no Gauss bound in the exponent applies.
DYADIC = ExponentialGrading(1.0, 2.0)


def graded_ends(a: float, b: float, toward_a: bool, toward_b: bool, layers: int, ratio: float) -> list[float]:
    r"""The interior panel ends of :math:`[a, b]`, ``layers`` geometric layers toward each graded end.

    A graded end at distance :math:`w\,\rho^j`, :math:`j = 1, \dots, L`, with
    :math:`w` the half width when both ends are graded and the whole width
    when one is. Shared by the panels toward walls and interfaces and the
    slab's grazing cosines (:mod:`.assembly`).
    """
    width = (b - a) / 2.0 if toward_a and toward_b else b - a
    depths = width * ratio ** np.arange(1, layers + 1)
    lower = list(a + depths[::-1]) if toward_a else []
    upper = list(b - depths) if toward_b else []
    return lower + upper


def exponential_ends(stop: np.ndarray, start: np.ndarray, steepest: np.ndarray, gentlest: np.ndarray, grading: ExponentialGrading) -> np.ndarray:
    r"""The ends of the intervals of ``[start, stop]`` graded exponentially toward ``stop``, sorted, ``(..., K + 2)``.

    The interval carries a family of attenuations :math:`e^{-\Sigma |{\rm stop} - s|}`
    with :math:`\Sigma` from ``gentlest`` to ``steepest`` (one member when the
    two are equal). The ends are the depths :math:`h\rho^k/\Sigma_{\rm steep}`,
    :math:`k = 0, \dots, K - 1`, from ``stop``, clipped to the interval, with
    :math:`h` and :math:`\rho` the ``grading``'s first depth and ratio
    (:class:`ExponentialGrading`),
    :math:`K = \lceil \log_\rho(D\,\Sigma_{\rm steep}/(h\,\Sigma_{\rm gentle})) \rceil + 1`
    and :math:`D` the :data:`VANISHING_DEPTH`: the first end is the steepest
    member's first piece, and the last reaches past :math:`D` mean free
    paths of the gentlest, so every member's own ends are resolved. None
    where the steepest member's optical width over the interval is at most
    :math:`h`.
    """
    width = np.abs(stop - start)
    thick = steepest * width > grading.first
    first = grading.first / np.where(thick, steepest, 1.0)
    spread = np.where(thick, steepest / np.where(thick, gentlest, 1.0), 1.0)
    steps = np.ceil(np.log(VANISHING_DEPTH * spread / grading.first) / np.log(grading.ratio)).astype(int) + 1   # (...,)
    k = np.arange(int(steps.max(initial=1)))
    depth = np.minimum(first[..., None] * grading.ratio ** k, width[..., None])
    depth = np.where(thick[..., None] & (k < steps[..., None]), depth, 0.0)
    inward = np.sign(start - stop)[..., None]
    return np.sort(np.concatenate([start[..., None], stop[..., None] + inward * depth, stop[..., None]], axis=-1), axis=-1)


#: The hp grading's tightness along a line: every piece no wider than its distance to the orbit coordinate's branch
#: points over it, so in the piece's coordinate on :math:`[-1, 1]` they lie at least at :math:`z = -1 \pm 2ic` from it
#: (the user's ruling of 2026-10-13, ``.claude/plans/characteristic_reference_architecture.md`` P3 step 3f: 9 to 13
#: nodes per piece at :math:`c = 2`, against up to 17 at :math:`c = 1`).
TIGHTNESS = 2.0

#: The error a piece's Gauss rule is held to, relative to the attenuation's integral over it: :func:`gauss_nodes`'s
#: target, 64 eps, the band the rule itself holds on the polynomials it is exact for (:mod:`.attenuation_rule`), so the
#: count spends its budget down to the rule's own rounding and no further (the user's ruling of 2026-10-13: a target of
#: 1e-14 sat below that floor, and PL1 read 1.09e-14 at it on the cylinder's speed family).
ALONG_LINE_ERROR = BAND

#: The Bernstein ellipses :func:`gauss_nodes` tries, as fractions of the way from 1 to the singularity's: each gives a
#: valid bound, and their least is still one.
_ELLIPSES = np.geomspace(0.02, 0.98, 24)


def branch_walk(length: np.ndarray, height: np.ndarray, approach: np.ndarray, tightness: float = TIGHTNESS) -> np.ndarray:
    r"""The hp grading of a slot toward its end nearer the closest approach: the piece ends' distances from that end, ``(..., K)``; 0 past the last.

    A slot of length ``length`` whose near end lies ``approach`` (:math:`e`,
    the parameter distance from the closest approach :math:`t^*`) from
    :math:`t^*`, on a line whose branch points sit at
    :math:`t^* \pm i\,` ``height`` (:math:`b/|P\Omega|`), is cut at
    :math:`a_{k+1} = a_k + D(a_k)/c`, :math:`a_0 = 0`, with
    :math:`D(a)^2 = (b/|P\Omega|)^2 + (e + a)^2` the distance from the point
    :math:`a` beyond the near end to the branch points. :math:`D` grows away
    from the closest approach, so it is each piece's least distance and every
    piece is no wider than it over :math:`c` (``tightness``). The distance is
    formed in closed form from exact quantities, never through the line's
    parameter, which rounds the end of a small cavity's slot onto the closest
    approach (`[M]` 2026-10-13, a cavity of radius 1e-12: the walk read
    :math:`D = 0` and stalled; ``scratch/characteristic_architecture/p3/product/qa4_fix.log``).

    Since :math:`D \ge e + a` and :math:`D \ge D(0) = d_0 > 0`, :math:`a_1 = d_0/c`
    and :math:`a_{k+1} \ge a_k (1 + 1/c)`, so the walk has at most
    :math:`1 + \lceil \ln(c(L + d_0)/d_0) / \ln(1 + 1/c) \rceil` steps. A slot
    with :math:`d_0 = 0` (a line through the centre, ending at it) is entire
    in the parameter and not cut.
    """
    length, height, approach = np.broadcast_arrays(*(np.asarray(a, dtype=float) for a in (length, height, approach)))
    nearest = np.hypot(height, approach)                                                   # d_0
    graded = (nearest > 0.0) & (length > 0.0)
    d0 = np.where(graded, nearest, 1.0)
    steps = np.where(graded, np.ceil(np.log(tightness * (length + d0) / d0) / np.log1p(1.0 / tightness)) + 1, 0).astype(int)
    a = np.zeros_like(length)
    ends = []
    for _ in range(int(steps.max(initial=0))):
        a = np.where(graded & (a < length), a + np.hypot(height, approach + a) / tightness, a)
        ends.append(np.where(graded & (a < length), a, 0.0))
    if np.any(graded & (a < length)):
        raise RuntimeError("the hp walk did not reach the slot's far end within its closed-form step bound")
    return np.stack(ends, axis=-1) if ends else np.zeros(length.shape + (0,))


def bernstein_parameter(z: np.ndarray) -> np.ndarray:
    r"""The parameter :math:`\rho = |z \pm \sqrt{z^2 - 1}|`, the larger, of the Bernstein ellipse through each complex ``z``: 1 on :math:`[-1, 1]`."""
    root = np.sqrt(np.asarray(z, dtype=complex) ** 2 - 1.0)
    return np.maximum(np.abs(z + root), np.abs(z - root))


def exact_nodes(degree: np.ndarray | int) -> np.ndarray:
    r"""The fewest nodes of a Gauss rule exact on every polynomial of degree ``degree``: :math:`\lceil (k + 1)/2 \rceil`."""
    return (np.asarray(degree, dtype=int) + 2) // 2


def gauss_nodes(singular: np.ndarray, majorant: Callable[[np.ndarray], np.ndarray], target: float = ALONG_LINE_ERROR) -> np.ndarray:
    r"""The least node count :math:`m` of a Gauss rule of a positive weight that integrates the factor to ``target``, per piece, ``(...)``.

    For a positive weight :math:`\omega` the :math:`m`-point Gauss rule's
    error on :math:`f\omega` is at most :math:`2\int\omega\;E_{2m-1}(f)`, with
    :math:`E_k(f)` the best uniform error of a polynomial of degree :math:`k`
    (the rule is exact on those and its weights are positive), and
    :math:`E_k(f) \le 2M(\rho)\rho^{-k}/(\rho - 1)` for any Bernstein ellipse
    :math:`\rho` on and inside which :math:`f` is analytic and bounded by
    :math:`M(\rho)` (Trefethen, *Approximation Theory and Approximation
    Practice*, 2013, Theorem 8.2). So :math:`m` is :func:`exact_nodes` of
    the least degree :math:`k` at which
    :math:`4M(\rho)\rho^{-k}/(\rho - 1) \le` ``target`` on one of 24
    ellipses between 1 and :math:`\rho_{\rm sing}` (:data:`_ELLIPSES`):
    each gives a valid bound, so the count is safe, though not the least
    the theorem allows. The count does not depend on the weight, so not on
    the optical width.

    ``singular`` is :math:`\rho_{\rm sing}`, the ellipse through each piece's
    nearest singularity (:func:`bernstein_parameter`), ``(...)``;
    ``majorant(rho)`` returns :math:`M(\rho)` for an ellipse ``rho`` ``(...)``.
    Every count is at least 1. Refused: a singularity on the interval
    (:math:`\rho_{\rm sing} \le 1`, no ellipse to bound on) and a majorant
    that is not finite and positive.
    """
    singular = np.asarray(singular, dtype=float)
    if not np.all(singular > 1.0):
        raise ValueError("a node count is bounded on a Bernstein ellipse inside the singularity's: it lies on the interval")
    degree = np.full(singular.shape, np.inf)
    for fraction in _ELLIPSES:
        rho = 1.0 + (singular - 1.0) * fraction
        bound = majorant(rho)
        if not np.all(np.isfinite(bound) & (bound > 0.0)):
            raise ValueError("a majorant of the integrand on an ellipse is finite and positive")
        needed = np.ceil(np.log(4.0 * bound / ((rho - 1.0) * target)) / np.log(rho))
        np.minimum(degree, np.maximum(needed, 1.0), out=degree)
    return exact_nodes(degree.astype(int))


def halvings(width: np.ndarray, distance: np.ndarray) -> np.ndarray:
    r"""Halvings of a piece of ``width`` until it is no wider than ``distance``, its distance to a singularity (hp).

    :math:`\lceil\log_2(w/d)\rceil`, 0 where the piece is already no wider,
    at most ``np.finfo(float).nmant``. Refused: a positive width at a distance
    of zero or less, a singularity on the interval itself.
    """
    width, distance = np.asarray(width, dtype=float), np.asarray(distance, dtype=float)
    needed = distance < width
    if np.any(needed & (distance <= 0.0)):
        raise ValueError("an hp grading toward a singularity on the interval itself has no finite depth")
    ratio = np.divide(width, distance, out=np.ones_like(width * distance), where=needed)
    return np.where(needed, np.minimum(np.ceil(np.log2(ratio)), np.finfo(float).nmant), 0).astype(int)


__all__ = [
    "ALONG_LINE_ERROR",
    "DYADIC",
    "ExponentialGrading",
    "TIGHTNESS",
    "VANISHING_DEPTH",
    "bernstein_parameter",
    "branch_walk",
    "exponential_ends",
    "exact_nodes",
    "gauss_nodes",
    "graded_ends",
    "halvings",
]
