r"""The gradings of the characteristic reference: where a quadrature puts its piece ends.

Three laws, each the one place it is spelled:

* **geometric** (:func:`graded_ends`): layers at :math:`w\rho^j` toward an
  end, the panels toward walls and interfaces and the halvings of the line
  rules;
* **exponential** (:func:`exponential_ends`): ends at a first depth times
  :math:`\rho^k` mean free paths toward where an attenuation concentrates,
  until it vanishes (:data:`VANISHING_DEPTH`), over a family of attenuations
  at once (the cylinder's polar angles, read as speeds along one in-plane
  line); the first depth and :math:`\rho` are derived from the node count
  of the rule integrating the pieces (:class:`ExponentialGrading`);
* **hp** (:func:`halvings`): halve toward a singularity off the interval
  until each piece is no wider than its distance to it, so Gauss-Legendre
  converges geometrically on every piece. The branch points of a line's
  orbit coordinate (the traversal rule), the next radius and :math:`b = 0`
  (the impact rule), and the grazing direction (the direction rules) are
  graded by it.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import e, lgamma, log, exp

import numpy as np

#: The optical depth past which an attenuation vanishes in double precision (e^-64): the last graded end reaches it.
VANISHING_DEPTH = 64.0

#: The widest first piece an exponential grading cuts, in mean free paths. The Gauss bound sees e^-x alone; the integrand
#: is a basis polynomial times it, whose variation over the panel the piece does not know, so the cap is measured, not
#: derived. `[M]` 2026-10-12: at a cap of 2 the cylinder's slowest polar node missed its mpmath density by 2.5e-13 at the
#: wall and 1.2e-13 at an interface at 8 points (the test-architect's P1 rows, ``scratch/characteristic_architecture/p3/
#: step3/ta/e/``), within its 1e-13 band at a cap of 1 (ruled 2026-10-12: the old law's first depth). Uncapped, the
#: nodes' own rounding (about h eps relative) also shows: 4 to 40 eps at 12 to 24 points, 7.9 to 32.6 mean free paths wide.
_WIDEST_FIRST = 1.0


@dataclass(frozen=True)
class ExponentialGrading:
    r"""Where an exponential grading puts its ends, in mean free paths of the attenuation: a first depth and a ratio.

    Each piece of the grading is integrated by an :math:`n`-point
    Gauss-Legendre rule, whose error on :math:`e^{-x}` over a piece
    :math:`[d, d + h]` is at most
    :math:`c_n h^{2n+1} e^{-d}`, :math:`c_n = (n!)^4 / ((2n+1)\,((2n)!)^3)`
    (Abramowitz and Stegun 25.4.30). Two widths follow from it, with
    :math:`h_n = (\epsilon / c_n)^{1/(2n+1)}` the widest piece at the
    attenuated end that :math:`n` points integrate to the double-precision
    :math:`\epsilon`:

    * the **first** piece :math:`[0, h]`, the attenuation's whole weight:
      :math:`h = \min(h_n, 1)` (:data:`_WIDEST_FIRST`, measured, since the
      integrand is a basis polynomial times :math:`e^{-x}`). A stretch no
      wider than it is not cut;
    * the **geometric** pieces :math:`[d, \rho d]`, whose bound
      :math:`c_n (d(\rho - 1))^{2n+1} e^{-d}` is largest at
      :math:`d = 2n + 1`: :math:`\rho = 1 + e\,h_n/(2n + 1)`.

    `[M]` 2026-10-12 (``scratch/characteristic_architecture/p3/step3/ratio/law_cap1.log``), the integral of
    :math:`e^{-x}` over :math:`[0, 64]` read at 0, in double precision:

    ====  =======  =======  ======  =======  ===========  ===================
    n     first    rho      pieces  error    ratio 2 [1]  ratio sqrt 2 [1]
    ====  =======  =======  ======  =======  ===========  ===================
    4     0.194    1.059    103     2.7e-15  1.1e-6       5.0e-9
    6     1        1.212    23      1.1e-16  9.7e-10      2.6e-13
    8     1        1.419    13      0        7.4e-13      2.2e-16
    12    1        1.856    8       0        0            1.1e-16
    16    1        2.238    7       0        0            1.1e-16
    ====  =======  =======  ======  =======  ===========  ===================

    [1] the fixed laws this replaced, first depth 1 (7 and 13 pieces): the
    error is the node count's, and no fixed ratio holds at every count
    (ruled 2026-10-12: derive it from the node count). Below 5 points the
    rounding over the many pieces is what remains (1.2e-13 at 2 points, 4522
    pieces).
    """

    first: float
    ratio: float

    @classmethod
    def of(cls, points: int) -> "ExponentialGrading":
        r"""The grading a rule of ``points`` Gauss-Legendre nodes per piece integrates to double precision."""
        if points < 1:
            raise ValueError(f"a Gauss-Legendre rule has at least one point; got {points}")
        n = points
        log_bound = 4.0 * lgamma(n + 1) - log(2 * n + 1) - 3.0 * lgamma(2 * n + 1)          # log c_n
        resolved = exp((log(np.finfo(float).eps) - log_bound) / (2 * n + 1))               # h_n
        return cls(min(resolved, _WIDEST_FIRST), 1.0 + e * resolved / (2 * n + 1))


#: The dyadic grading, first depth 1 and ratio 2: for a coordinate in which the attenuation's exponent is not linear (the
#: polar angle and the cosine toward the normal, :mod:`.lines`), where the Gauss bound in the exponent does not apply.
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
    (:class:`ExponentialGrading`, from the integrating rule's node count),
    :math:`K = \lceil \log_\rho(D\,\Sigma_{\rm steep}/(h\,\Sigma_{\rm gentle})) \rceil + 1`
    and :math:`D` the :data:`VANISHING_DEPTH`: the first end is the steepest
    member's first piece, and the last reaches past :math:`D` mean free
    paths of the gentlest, so every member's own ends are resolved. None
    where the steepest member's optical width over the interval is at most
    :math:`h`.

    `[M]` 2026-10-11 (the test-architect, P3 step 3a, ``scratch/characteristic_architecture/p3/step3/ta/``):
    on the cylinder's polar family a first end 16 times the steepest mean free
    path cost 6.7e-11 on a thick wall at 8 points and 100 times cost 1.3e-6;
    stopping at 7 doublings of the steepest missed every closed-form reading
    by 7.8e-9 to 2.7e-3.
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


__all__ = ["DYADIC", "ExponentialGrading", "VANISHING_DEPTH", "exponential_ends", "graded_ends", "halvings"]
