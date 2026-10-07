r"""The panel basis: the space the emission density and the flux are represented in.

Each region :math:`[r_k, r_{k+1}]` of the body is cut into **panels**, graded
geometrically toward each end that is a wall or an interface, where the flux
has its boundary layer and its derivative singularities. A singular stratum
(a solid body's centre or axis) is not graded: the flux is even and smooth
there. On each panel the basis is the :math:`p + 1` nodal Lagrange functions
through the panel's Gauss-Legendre points, so it is discontinuous across
panel ends and no node lies on one; a coefficient is a value at a node, which
is what lets the dense pencil's single-sign test read coefficients as values.

**Even at a singular stratum.** On the panel whose lower end is the centre
or the axis, the Lagrange functions are polynomials in :math:`c^2` through
the squares of the same nodes: they span :math:`1, c^2, \dots, c^{2p}`.
A smooth function invariant under the stratum's isotropy :math:`O(d)` is a
smooth function of :math:`|x|^2` (Schwarz's theorem), so the flux there has
no odd mode, and an odd mode is what the basis would add: the line integral
of :math:`c^{2m+1}` is an Abel transform carrying a :math:`b^{2m+2}\log b`
term, which held the closed sphere's centre panel to algebraic convergence
in the impact parameter (measured 2026-10-06: 3.5e-8, 1.6e-10, 7.0e-13 at 8, 16,
32 points; the user's ruling the same day).

**The panel ends are a partition.** They are posed as a
:class:`~orpheus.geometry.chord.ConcentricPartition`, a refinement of the
body's with the same ends, so the geometric kernel's chord through it yields
every piece of every line, one slot per panel crossed, with its
cancellation-free length and its panel read from the slot's region code (the
user's ruling of 2026-10-06, P1 step (b) second rung).

**The mass matrix** is :math:`W_{ij} = \int u_i u_j\,\mathrm{d}V` in the
chart's volume measure, whose density in the orbit coordinate,
:math:`c\,d\,r^{d-1}`, is the kernel's derivative of the one definition of
the measure (:meth:`~orpheus.geometry.chart.Chart.measure_density`). The
integrand is a polynomial of degree at most :math:`4p + d - 1` (on the even
panel), so Gauss-Legendre with :math:`2p + 2` points on every panel
integrates every entry exactly.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from typing import NamedTuple

import numpy as np
from scipy.linalg import block_diag

from orpheus.derivations.common.quadrature import composite_gauss_legendre, gauss_legendre
from orpheus.geometry.chord import ConcentricPartition

from .grading import graded_ends


def _on_stratum(regions: ConcentricPartition, orbit_coordinate: np.ndarray) -> np.ndarray:
    """Whether each orbit coordinate is a singular stratum of the body's chart (a solid body's centre or axis)."""
    return np.isin(orbit_coordinate, [s.orbit_value for s in regions.chart.singular_strata])


def _even_coordinate(reference: np.ndarray) -> np.ndarray:
    r"""An even panel's coordinate :math:`(c/h)^2` from the reference :math:`x = 2c/h - 1`: :math:`((x + 1)/2)^2`.

    The one map for the points and the nodes, so the basis interpolates at its own nodes.
    """
    return ((reference + 1.0) / 2.0) ** 2


class _LagrangeTables(NamedTuple):
    r"""The two panel coordinates' Lagrange tables, indexed by the coordinate: 0 the reference, 1 the even one.

    Attributes
    ----------
    nodes:
        The nodes :math:`\xi_k`, ``(2, p + 1)``.
    others:
        For each function :math:`i`, the indices :math:`k \ne i`, ``(p + 1, p)``.
    spread:
        :math:`\xi_i - \xi_k` over those :math:`k`, ``(2, p + 1, p)``.
    """

    nodes: np.ndarray
    others: np.ndarray
    spread: np.ndarray


@dataclass(frozen=True, eq=False)
class PanelBasis:
    r"""Discontinuous nodal Lagrange panels of degree :math:`p` on a concentric body.

    Built by :meth:`of`. Node :math:`i` is the :math:`m`-th Gauss-Legendre point
    of panel :math:`P`, :math:`i = P(p + 1) + m`.

    Attributes
    ----------
    regions:
        The body's partition.
    partition:
        The panel ends, a refinement of ``regions`` with the same ends.
    degree:
        The polynomial degree :math:`p` on each panel.
    """

    regions: ConcentricPartition
    partition: ConcentricPartition
    degree: int

    def __post_init__(self) -> None:
        coarse, fine = np.asarray(self.regions.breakpoints), np.asarray(self.partition.breakpoints)
        if self.partition.chart != self.regions.chart:
            raise ValueError("the panel partition and the body's partition are posed on one chart")
        if fine[0] != coarse[0] or fine[-1] != coarse[-1] or not np.isin(coarse, fine).all():
            raise ValueError(
                f"the panel ends refine the body's breakpoints {tuple(coarse)} with the same ends; got {tuple(fine)}"
            )
        if self.degree < 0:
            raise ValueError(f"the panel degree is non-negative; got {self.degree}")

    @classmethod
    def of(cls, regions: ConcentricPartition, degree: int, layers: int, ratio: float) -> "PanelBasis":
        r"""Panels of degree ``degree`` graded by ``layers`` layers of ratio ``ratio`` toward every wall and interface."""
        if layers < 0 or not 0.0 < ratio < 1.0:
            raise ValueError(f"grading takes layers >= 0 and a ratio in (0, 1); got {layers}, {ratio}")
        r = np.asarray(regions.breakpoints)
        graded = ~_on_stratum(regions, r)
        ends = [r[0]]
        for k in range(len(r) - 1):
            ends += graded_ends(r[k], r[k + 1], graded[k], graded[k + 1], layers, ratio) + [r[k + 1]]
        partition = ConcentricPartition(regions.chart, tuple(ends), regions.pose)
        return cls(regions, partition, degree)

    # ── the nodes ────────────────────────────────────────────────────────

    @cached_property
    def region_of_panel(self) -> np.ndarray:
        """The region holding each panel, ``(P,)``."""
        coarse, fine = np.asarray(self.regions.breakpoints), np.asarray(self.partition.breakpoints)
        return np.searchsorted(coarse, fine[:-1], side="right") - 1

    @property
    def n_panels(self) -> int:
        """The number :math:`P` of panels."""
        return self.partition.n_regions

    @property
    def per_panel(self) -> int:
        """The number :math:`p + 1` of functions on each panel."""
        return self.degree + 1

    @property
    def size(self) -> int:
        """The dimension :math:`N = P(p + 1)` of the basis."""
        return self.n_panels * self.per_panel

    @cached_property
    def _reference_nodes(self) -> np.ndarray:
        return gauss_legendre(-1.0, 1.0, self.per_panel).pts

    @cached_property
    def even(self) -> np.ndarray:
        r"""Whether each panel's functions are polynomials in :math:`c^2`, ``(P,)``: the panel whose lower end is a singular stratum."""
        return _on_stratum(self.regions, np.asarray(self.partition.breakpoints[:-1]))

    @cached_property
    def nodes(self) -> np.ndarray:
        """The orbit coordinate of each node, ``(N,)``, in panel order."""
        return composite_gauss_legendre(self.partition.breakpoints, self.per_panel).pts

    @property
    def panel(self) -> np.ndarray:
        """The panel of each node, ``(N,)``."""
        return np.repeat(np.arange(self.n_panels), self.per_panel)

    @property
    def region(self) -> np.ndarray:
        """The region of each node, ``(N,)``."""
        return self.region_of_panel[self.panel]

    def on_panels(self, per_region: np.ndarray) -> np.ndarray:
        """A per-region table ``(n, ...)`` read on each panel, ``(P, ...)``."""
        per_region = np.asarray(per_region)
        if per_region.shape[:1] != (self.regions.n_regions,):
            raise ValueError(
                f"a per-region table has one entry per region, leading shape ({self.regions.n_regions},); "
                f"got {per_region.shape}"
            )
        return per_region[self.region_of_panel]

    # ── the functions ────────────────────────────────────────────────────

    def values(self, orbit_coordinate: np.ndarray, panel: np.ndarray) -> np.ndarray:
        r"""The :math:`p + 1` functions of ``panel`` at ``orbit_coordinate``, ``(..., p + 1)``.

        The Lagrange product form in the panel's coordinate, which is exact
        at a node (one-hot) and needs no division by a node distance. The
        coordinate is the reference :math:`x \in [-1, 1]`, or on an
        :attr:`even` panel :math:`[0, h]` its square :math:`(c/h)^2`, with the
        nodes squared alike.
        """
        ends = np.asarray(self.partition.breakpoints)
        panel = np.asarray(panel)
        a, b = ends[panel], ends[panel + 1]
        x = (2.0 * np.asarray(orbit_coordinate, dtype=float) - (a + b)) / (b - a)
        even = self.even[panel]
        x = np.where(even, _even_coordinate(x), x)
        nodes, others, spread = self._lagrange
        kind = even.astype(int)
        factors = (x[..., None] - nodes[kind])[..., others] / spread[kind]                 # (..., p + 1, p)
        return factors.prod(axis=-1)

    @cached_property
    def _lagrange(self) -> _LagrangeTables:
        r"""The Lagrange tables of the reference and the :attr:`even` coordinate. Built once, they leave :meth:`values` one product per point, which had been
        2.4 times slower rebuilding them at every point (`[M]` 2026-10-07,
        #586: 4e5 points of a three-region cylinder's basis).
        """
        reference = self._reference_nodes
        nodes = np.stack([reference, _even_coordinate(reference)])
        n = self.per_panel
        others = np.array([[k for k in range(n) if k != i] for i in range(n)], dtype=int).reshape(n, n - 1)
        spread = nodes[:, :, None] - nodes[:, others]
        return _LagrangeTables(nodes, others, spread)

    def columns(self, panel: np.ndarray) -> np.ndarray:
        """The basis indices of the functions of ``panel``, ``(..., p + 1)``."""
        return np.asarray(panel)[..., None] * self.per_panel + np.arange(self.per_panel)

    # ── the metric ───────────────────────────────────────────────────────

    @cached_property
    def mass(self) -> np.ndarray:
        r"""The Gram matrix :math:`W_{ij} = \int u_i u_j\,\mathrm{d}V`, ``(N, N)``, block-diagonal by panel."""
        points = 2 * self.per_panel
        rule = composite_gauss_legendre(self.partition.breakpoints, points)
        panel = np.repeat(np.arange(self.n_panels), points)
        shape = (self.n_panels, points)
        u = self.values(rule.pts, panel).reshape(*shape, self.per_panel)
        weight = (rule.wts * self.regions.chart.measure_density(rule.pts)).reshape(shape)
        return block_diag(*np.einsum("Pqi,Pq,Pqj->Pij", u, weight, u))


__all__ = ["PanelBasis"]
