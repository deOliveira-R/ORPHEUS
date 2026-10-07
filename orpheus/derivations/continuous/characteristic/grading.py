r"""The gradings of the characteristic reference: where a quadrature puts its piece ends.

Three laws, each the one place it is spelled:

* **geometric** (:func:`graded_ends`): layers at :math:`w\rho^j` toward an
  end, the panels toward walls and interfaces and the halvings of the line
  rules;
* **exponential** (:func:`exponential_ends`): ends at :math:`2^k` mean free
  paths toward where an attenuation concentrates, until it vanishes
  (:data:`VANISHING_DEPTH`);
* **hp** (:func:`halvings`): halve toward a singularity off the interval
  until each piece is no wider than its distance to it, so Gauss-Legendre
  converges geometrically on every piece. The branch points of a line's
  orbit coordinate (the traversal rule), the next radius and :math:`b = 0`
  (the impact rule), and the grazing direction (the direction rules) are
  graded by it.
"""

from __future__ import annotations

import numpy as np

#: Piece ends at 2^k mean free paths, k = 0..6, from a graded end: beyond 64, e^-64 is below double precision.
_DOUBLINGS = 7

#: The optical depth past which an attenuation vanishes in double precision: the last graded end, 2^6 = 64.
VANISHING_DEPTH = 2.0 ** (_DOUBLINGS - 1)

#: An interval of optical width at most this is one Gauss panel: e^-2 is integrated to machine precision.
_THIN = 2.0


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


def exponential_ends(stop: np.ndarray, start: np.ndarray, sigma: np.ndarray) -> np.ndarray:
    r"""The ends of the intervals of ``[start, stop]`` graded exponentially toward ``stop``, sorted, ``(..., K + 2)``.

    Depths :math:`2^k/\Sigma` from ``stop``, clipped to the interval; none where
    the interval's optical width is at most :data:`_THIN`.
    """
    width = np.abs(stop - start)
    thick = sigma * width > _THIN
    mean_free_path = 1.0 / np.where(thick, sigma, 1.0)
    depth = np.where(thick[..., None], np.minimum(mean_free_path[..., None] * 2.0 ** np.arange(_DOUBLINGS), width[..., None]), 0.0)
    inward = np.sign(start - stop)[..., None]
    return np.sort(np.concatenate([start[..., None], stop[..., None] + inward * depth, stop[..., None]], axis=-1), axis=-1)


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


__all__ = ["VANISHING_DEPTH", "exponential_ends", "graded_ends", "halvings"]
