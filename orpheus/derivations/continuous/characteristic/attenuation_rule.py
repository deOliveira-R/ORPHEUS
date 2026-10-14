r"""The Gauss rule of the attenuation weight: the quadrature along a piece of a line, with the exponential as its measure.

Along a piece of a line, every integral the transport forms is a smooth
factor (a panel function of the orbit coordinate, or a product of two)
against the attenuation to one of the piece's ends. In the piece's
reference coordinate :math:`s \in [-1, 1]`, with :math:`s = 1` the
attenuated end and :math:`x` the piece's optical half-width, that
attenuation is the weight

.. math::

    \omega_x(s) = e^{-x(1 - s)}, \qquad \int_{-1}^{1} \omega_x\,\mathrm ds = \frac{1 - e^{-2x}}{x} = 2\,\tilde\imath_0(x),

with :math:`\tilde\imath_0(x) = e^{-x} i_0(x)` the scaled modified spherical
Bessel function of order 0. :class:`AttenuationRule` is the :math:`m`-point
Gauss rule of :math:`\omega_x`: :math:`m` nodes and :math:`m` positive
weights that integrate :math:`p\,\omega_x` exactly for every polynomial
:math:`p` of degree at most :math:`2m - 1`, at every :math:`x \ge 0`. The
smooth factor is evaluated at the rule's own nodes, so nothing is
interpolated, and the attenuation costs no piece mesh however thick the
piece. Positive weights round componentwise: the rule's error on a
non-negative factor is a relative one.

**The nodes are offsets.** Each node is held as its offset
:math:`\delta = 1 - s \in (0, 2)` from the attenuated end, never as
:math:`s`. On a thick piece the nodes crowd at :math:`s = 1`, where
:math:`1 - s` rounded from a double :math:`s` keeps only
:math:`\epsilon / \delta` relative digits; the consumer places a node a
distance :math:`\delta h/2` from the piece's end (:meth:`~orpheus.geometry.chord.Rays.advanced`),
so the offset is the quantity it needs. `[M]` 2026-10-13
(``scratch/characteristic_architecture/p3/product/rb3_compare.log``):
the same rule held in :math:`s` integrates the Legendre polynomials to 93
eps of the zeroth moment at :math:`m = 20, x = 1000`, against 42 eps held
as offsets.

**The construction** (the discretised Stieltjes procedure, Gautschi,
*Orthogonal Polynomials: Computation and Approximation*, 2004, §2.2.3,
then Golub and Welsch, 1969). The weight is first discretised, a
quadrature of the *measure*, exact to rounding on the polynomials of
degree :math:`2m - 1` times it; the Stieltjes procedure then builds the
Jacobi matrix of that discrete measure, with normalised vectors; its eigenvalues are the nodes and
the squared first components of its eigenvectors the weights.

* The measure is written in :math:`v = \lambda\delta` (``_scale``):
  :math:`\lambda = x` where :math:`2x >` :data:`PANEL_WIDTH`, so the weight
  is :math:`e^{-v}\,\mathrm dv/x`, and :math:`\lambda = 1` below, so the
  variable is the offset itself, with no division by :math:`x`
  (:math:`x = 0` and a subnormal :math:`x` are exact).
* It is truncated at :data:`TRUNCATION` and discretised by panels of width
  at most :data:`PANEL_WIDTH`, :math:`m +` :data:`EXTRA_POINTS` points each.
  The recurrence runs on :math:`v`, so the nodes are offsets from the
  attenuated end and do not crowd.

The weights are rescaled to the exact zeroth moment, :math:`-\mathrm{expm1}(-2x)/x`,
so a constant is integrated to its last ulp.

There is no switch between constructions, and no table: the rule of the
truncated measure past :math:`2x \ge` :data:`TRUNCATION` is built the same
way. `[M]` 2026-10-13 (``.../p3/product/table_vs_stieltjes.log``): a table
of those rules generated in mpmath agreed with this construction to noise
(within 18 eps either way), so it was a second definition of one quantity
and was not shipped.

`[M]` 2026-10-13 (``.../p3/product/rb_batched_noconst.log``, this
construction, 45 to 48 values of :math:`x` per :math:`m`, from 0 and
5e-324 to 1000 with 31 in :math:`[m/2, 2m]`): every Legendre polynomial of
degree at most :math:`2m - 1` integrated to within 47.9 eps of the zeroth
moment for :math:`m = 1` to 23.

The design and its rulings: ``.claude/plans/characteristic_reference_architecture.md``
(P3 step 3f).
"""

from __future__ import annotations

from dataclasses import dataclass
from math import ceil, log

import numpy as np

from orpheus.derivations.common.quadrature import gauss_legendre

#: The offset :math:`u = x\delta` past which the weight's measure is dropped: :math:`\ln(1/\epsilon) + 10\ln 2`. Every
#: Legendre polynomial is bounded by 1 on :math:`[-1, 1]`, so the dropped tail moves each moment of the rule by at most
#: :math:`e^{-U}` relative to the zeroth, :math:`2^{-10}\epsilon`.
TRUNCATION = log(1.0 / np.finfo(float).eps) + 10.0 * log(2.0)

#: The widest panel of the measure's discretisation, in the offset :math:`u`; also the width :math:`2x` below which the
#: rule is built on one panel in the offset :math:`\delta` itself. Over such a panel :math:`m +` :data:`EXTRA_POINTS`
#: Gauss-Legendre points leave the polynomials of degree :math:`2m - 1` times the exponential 20 degrees of exactness
#: for the exponential's own variation, :math:`e^{4}` at most across a panel. The pair is measured, not derived: `[M]`
#: 2026-10-13 (``scratch/characteristic_architecture/p3/product/gate_spec.md`` R2, the arms ``wide panel`` and ``too few
#: extra points`` of the test-architect's battery) a wider panel or fewer points reddens RB3.
PANEL_WIDTH = 4.0

#: The Gauss-Legendre points of a discretisation panel beyond the rule's :math:`m` (:data:`PANEL_WIDTH`).
EXTRA_POINTS = 10

#: The band the rule holds on every polynomial it is exact for, relative to its zeroth moment, is 64 eps (RB3, `[M]`
#: 2026-10-13, ``.../p3/product/table_vs_stieltjes.log``, worst 61.4 eps over m = 1 to 40). It is the node law's target
#: (:data:`~.grading.ALONG_LINE_ERROR`).
BAND = 64.0 * float(np.finfo(float).eps)

#: The most nodes a rule is built with is the edge of the measured band. `[M]` 2026-10-13 (the numerics-investigator,
#: ``.../p3/product/table_vs_stieltjes.log``): within 64 eps of the zeroth moment on every polynomial it is exact for, at
#: every m from 1 to 40 (worst 61.4 eps, m = 39); a sparse check at m = 64 read 68.5 eps (qa's review).
#:
#: SCOPE-BOUNDARY[guard] machinery: a construction held inside the band past 40 nodes (a wider discretisation, or the
#: rule's error measured and carried as its own band). ruling: the user, 2026-10-13, "Refuse above the measured m",
#: ``.claude/plans/characteristic_reference_architecture.md`` P3 step 3f. revisit: a basis degree whose node law asks
#: more than 40 (degree 16 asks 41 on the triangle; the working degree 3 asks 10 to 11).
MOST_NODES = 40

#: The most floats one batch of rules holds in its discretised measures, ``(batch, discretisation)``: the batch is cut to it.
_MEASURE_BUDGET = 1 << 24


def zeroth_moment(x: np.ndarray) -> np.ndarray:
    r""":math:`\int_{-1}^{1} e^{-x(1 - s)}\,\mathrm ds = -\mathrm{expm1}(-2x)/x`, 2 at :math:`x = 0`; ``x`` non-negative."""
    x = np.asarray(x, dtype=float)
    positive = x > 0.0
    return np.where(positive, -np.expm1(-2.0 * x) / np.where(positive, x, 1.0), 2.0)


@dataclass(frozen=True, eq=False)
class TriangleRule:
    r"""The rule of the weight :math:`e^{-x(s - s')}` on the triangle :math:`-1 \le s' \le s \le 1`, a lead and a lag point per node.

    :math:`\int_{-1}^{1} f(s) \int_{-1}^{s} g(s')\,e^{-x(s - s')}\,\mathrm ds'\,\mathrm ds
    \approx \sum_n w_n f(s_n) g(s'_n)`: the source :math:`g` at the lag point
    :math:`s'`, upstream, attenuated to the lead point :math:`s`, where it
    meets :math:`f`. Built by :meth:`AttenuationRule.triangle`. Both points
    are held as offsets from the piece's attenuated end :math:`s = 1`, like
    the one-level rule's nodes.

    Attributes
    ----------
    lead:
        The offset :math:`1 - s` of each node's lead point, ``(..., m, m)``.
    lag:
        The offset :math:`1 - s'` of each node's lag point, ``(..., m, m)``; the lead's plus the lag
        :math:`\tau = s - s'`.
    weight:
        The node's weight, ``(..., m, m)``, non-negative.
    """

    lead: np.ndarray
    lag: np.ndarray
    weight: np.ndarray


@dataclass(frozen=True, eq=False)
class AttenuationRule:
    r"""The Gauss rule of :math:`e^{-x(1 - s)}` on :math:`[-1, 1]`, per half-width :math:`x`, nodes held as offsets from :math:`s = 1`.

    :math:`\int_{-1}^{1} f(s)\,e^{-x(1 - s)}\,\mathrm ds \approx \sum_q w_q f(1 - \delta_q)`,
    exact for :math:`f` a polynomial of degree at most :math:`2m - 1`. Toward
    the start :math:`s = -1`, the weight :math:`e^{-x(1 + s)}`, the same rule
    is read from the other end: :math:`f(-1 + \delta_q)`.

    Attributes
    ----------
    offset:
        :math:`\delta_q = 1 - s_q`, ``(..., m)``, in :math:`(0, 2)`; on padding, the rule's last offset.
    weight:
        :math:`w_q`, ``(..., m)``, positive; 0 on padding.
    count:
        The live nodes of each rule, ``(...)``: the first ``count`` of the ``m``.
    """

    offset: np.ndarray
    weight: np.ndarray
    count: np.ndarray

    @classmethod
    def of(cls, half_width: np.ndarray, count: np.ndarray | int) -> "AttenuationRule":
        r"""The rule of ``count`` nodes at each optical half-width ``half_width`` (non-negative), padded to the largest count.

        ``count`` broadcasts against ``half_width``, ``(...)``, each at least 1.
        The rules of one count are built together.
        """
        x = np.asarray(half_width, dtype=float)
        counts = np.broadcast_to(np.asarray(count, dtype=int), x.shape)
        if not np.all(np.isfinite(x)) or np.any(x < 0.0):
            raise ValueError("an attenuation's optical half-width is finite and non-negative")
        if np.any(counts < 1):
            raise ValueError(f"a Gauss rule has at least one node; got {int(counts.min())}")
        if np.any(counts > MOST_NODES):
            raise ValueError(
                f"an attenuation rule is measured within its band up to {MOST_NODES} nodes; {int(counts.max())} asked "
                "(a scope boundary: MOST_NODES)"
            )
        width = int(counts.max(initial=1))
        offset = np.zeros(x.shape + (width,))
        weight = np.zeros(x.shape + (width,))
        for m in np.unique(counts):
            chosen = counts == m
            delta, w = _built(x[chosen], int(m))
            offset[chosen] = np.concatenate([delta, np.repeat(delta[:, -1:], width - m, axis=-1)], axis=-1)
            weight[chosen, :m] = w
        return cls(offset, weight, counts.copy())

    def triangle(self) -> TriangleRule:
        r"""The rule of :math:`e^{-x(s - s')}` on the triangle :math:`s' \le s` (:class:`TriangleRule`), from this rule.

        With the lag :math:`\tau = s - s' \in [0, 2]` the triangle is
        :math:`\int_0^2 e^{-x\tau} G(\tau)\,\mathrm d\tau`,
        :math:`G(\tau) = \int_{-1}^{1 - \tau} f(s' + \tau)\,g(s')\,\mathrm ds'`:
        the outer integral is this rule, whose offsets are the lags
        :math:`\tau_q = \delta_q` (the weight :math:`e^{-x(1 - s)}` read at
        :math:`\tau = 1 - s`); the inner one has no exponential and is
        Gauss-Legendre of the rule's width over :math:`[-1, 1 - \tau_q]`. The
        lead point's offset is :math:`(2 - \tau_q)(1 - g_k)/2`, formed with
        no difference near the end, and the lag's is the lead's plus
        :math:`\tau_q`.
        """
        m = self.offset.shape[-1]
        legendre = gauss_legendre(-1.0, 1.0, m)
        half = (2.0 - self.offset)[..., None] / 2.0                                         # (..., m, 1)
        lead = half * (1.0 - legendre.pts)                                                   # (..., m, m)
        return TriangleRule(lead, lead + self.offset[..., None], self.weight[..., None] * half * legendre.wts)


def _scale(x: np.ndarray) -> np.ndarray:
    r"""The scale :math:`\lambda` of the variable :math:`v = \lambda\delta` each rule is built in: :math:`x` where :math:`2x >` :data:`PANEL_WIDTH`, else 1.

    In :math:`v = x\delta` the weight is :math:`e^{-v}`, whose discretisation
    no longer depends on :math:`x` and whose nodes do not crowd; below the
    panel width the offset itself, with no division by :math:`x`, so
    :math:`x = 0` and a subnormal :math:`x` are exact.
    """
    return np.where(2.0 * x > PANEL_WIDTH, x, 1.0)


def _built(x: np.ndarray, m: int) -> tuple[np.ndarray, np.ndarray]:
    r"""The ``m``-node rules at half-widths ``x`` ``(B,)``: offsets and weights ``(B, m)``, weights at the exact zeroth moment.

    The rules of one discretisation (one panel count) are built together, and
    each batch is cut to :data:`_MEASURE_BUDGET`.
    """
    offset, weight = np.empty(x.shape + (m,)), np.empty(x.shape + (m,))
    panels = np.ceil(np.minimum(2.0 * _scale(x), TRUNCATION) / PANEL_WIDTH).astype(int)
    for k in np.unique(panels):
        chosen = np.flatnonzero(panels == k)
        size = (m + EXTRA_POINTS) * int(k)
        for batch in np.array_split(chosen, max(1, ceil(chosen.size * size / _MEASURE_BUDGET))):
            offset[batch], weight[batch] = _stieltjes(x[batch], m, int(k))
    weight *= (zeroth_moment(x) / weight.sum(axis=-1))[:, None]
    return offset, weight


def _stieltjes(x: np.ndarray, m: int, panels: int) -> tuple[np.ndarray, np.ndarray]:
    r"""The rules at half-widths ``x`` ``(B,)`` sharing one discretisation of ``panels`` panels.

    The measure :math:`e^{-x\delta}\,\mathrm d\delta = e^{-(x/\lambda)v}\,\mathrm dv/\lambda`
    in :math:`v = \lambda\delta` (``_scale``), on :math:`[0, \min(2\lambda,`
    :data:`TRUNCATION` :math:`)]`. Returns the offsets and the Christoffel
    weights, ``(B, m)``, before the rescale.
    """
    legendre = gauss_legendre(-1.0, 1.0, m + EXTRA_POINTS)
    scale = _scale(x)[:, None]
    top = np.minimum(2.0 * scale, TRUNCATION)
    ends = np.linspace(0.0, 1.0, panels + 1) * top                                            # (B, panels + 1)
    lower, upper = ends[:, :-1, None], ends[:, 1:, None]
    support = (lower + (upper - lower) * (legendre.pts + 1.0) / 2.0).reshape(x.size, -1)
    mass = ((upper - lower) / 2.0 * legendre.wts).reshape(x.size, -1) * np.exp(-(x[:, None] / scale) * support) / scale
    diagonal, off_diagonal = _jacobi(support, mass, m)
    jacobi = np.zeros((x.size, m, m))
    i = np.arange(m)
    jacobi[:, i, i] = diagonal
    jacobi[:, i[:-1], i[1:]] = jacobi[:, i[1:], i[:-1]] = off_diagonal
    nodes, vectors = np.linalg.eigh(jacobi)
    christoffel = mass.sum(axis=-1)[:, None] * vectors[:, 0, :] ** 2
    return nodes / scale, christoffel


def _jacobi(support: np.ndarray, mass: np.ndarray, m: int) -> tuple[np.ndarray, np.ndarray]:
    r"""The Jacobi matrix of each discrete measure (``support``, ``mass``) ``(B, N)``: its diagonal ``(B, m)`` and off-diagonal ``(B, m - 1)``.

    The Stieltjes procedure on normalised vectors (the Lanczos form), with no
    re-orthogonalisation: `[M]` 2026-10-13 (the test-architect's battery,
    ``scratch/characteristic_architecture/p3/product/ta/carve/``), a full
    re-orthogonalisation pass per step reddened no gate up to :math:`m = 80`,
    past :data:`MOST_NODES` (the user's ruling the same day: removed).
    """
    batch = support.shape[0]
    diagonal, off_diagonal = np.zeros((batch, m)), np.zeros((batch, m - 1))
    current = np.ones_like(support) / np.sqrt(mass.sum(axis=-1))[:, None]
    previous, coupling = np.zeros_like(support), np.zeros(batch)
    for k in range(m):
        diagonal[:, k] = np.sum(mass * support * current * current, axis=-1)
        if k == m - 1:
            break
        residual = (support - diagonal[:, k, None]) * current - coupling[:, None] * previous
        coupling = np.sqrt(np.sum(mass * residual * residual, axis=-1))
        off_diagonal[:, k] = coupling
        previous, current = current, residual / coupling[:, None]
    return diagonal, off_diagonal


__all__ = ["BAND", "EXTRA_POINTS", "MOST_NODES", "PANEL_WIDTH", "TRUNCATION", "AttenuationRule", "TriangleRule", "zeroth_moment"]
