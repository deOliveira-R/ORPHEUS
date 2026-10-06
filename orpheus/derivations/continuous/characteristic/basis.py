r"""The panel basis: the space the emission density and the flux are represented in.

Each region :math:`[r_k, r_{k+1}]` of the body is cut into **panels**, graded
geometrically toward each end that is a wall or an interface, where the flux
has its boundary layer and its derivative singularities. A singular stratum
(a solid body's centre or axis) is not graded: the flux is even and smooth
there. On each panel the basis is the :math:`p + 1` nodal Lagrange functions
through the panel's Gauss-Legendre points, so it is discontinuous across
panel ends and no node lies on one; a coefficient is a value at a node, which
is what lets the dense pencil's single-sign test read coefficients as values.

**The panel ends are a partition.** They are posed as a
:class:`~orpheus.geometry.chord.ConcentricPartition`, a refinement of the
body's with the same ends, so the geometric kernel's chord through it yields
every piece of every line, one slot per panel crossed, with its
cancellation-free length and its panel read from the slot's region code (the
user's ruling of 2026-10-06, P1 step (b) second rung).

**The mass matrix** is :math:`W_{ij} = \int u_i u_j\,\mathrm{d}V` in the
chart's volume measure, whose density in the orbit coordinate is derived
here from the one definition of the measure,
:math:`m = c\,(T(r_{j+1}) - T(r_j))` with :math:`T(r) = r^d`
(:meth:`~orpheus.geometry.coord.CoordSystem.measure`): the density is
:math:`c\,d\,r^{d-1}`, a polynomial of degree :math:`d - 1 \le 2`, so
Gauss-Legendre with :math:`p + 2` points per panel integrates every entry
exactly.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

import numpy as np
from scipy.linalg import block_diag

from orpheus.derivations.common.quadrature import composite_gauss_legendre, gauss_legendre
from orpheus.geometry.chord import ConcentricPartition


def _graded_ends(a: float, b: float, toward_a: bool, toward_b: bool, layers: int, ratio: float) -> list[float]:
    r"""The interior panel ends of :math:`[a, b]`, ``layers`` geometric layers toward each graded end.

    A graded end at distance :math:`w\,\rho^j`, :math:`j = 1, \dots, L`, with
    :math:`w` the half width when both ends are graded and the whole width
    when one is.
    """
    width = (b - a) / 2.0 if toward_a and toward_b else b - a
    depths = width * ratio ** np.arange(1, layers + 1)
    lower = list(a + depths[::-1]) if toward_a else []
    upper = list(b - depths) if toward_b else []
    return lower + upper


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
        strata = {s.orbit_value for s in regions.chart.singular_strata}
        graded = [value not in strata for value in r]
        ends = [r[0]]
        for k in range(len(r) - 1):
            ends += _graded_ends(r[k], r[k + 1], graded[k], graded[k + 1], layers, ratio) + [r[k + 1]]
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

        The Lagrange product form, which is exact at a node (one-hot) and
        needs no division by a node distance.
        """
        ends = np.asarray(self.partition.breakpoints)
        panel = np.asarray(panel)
        a, b = ends[panel], ends[panel + 1]
        x = (2.0 * np.asarray(orbit_coordinate, dtype=float) - (a + b)) / (b - a)
        xi = self._reference_nodes
        others = ~np.eye(self.per_panel, dtype=bool)
        spread = np.where(others, xi[:, None] - xi[None, :], 1.0)
        factors = np.where(others, (x[..., None, None] - xi[None, :]) / spread, 1.0)
        return factors.prod(axis=-1)

    def columns(self, panel: np.ndarray) -> np.ndarray:
        """The basis indices of the functions of ``panel``, ``(..., p + 1)``."""
        return np.asarray(panel)[..., None] * self.per_panel + np.arange(self.per_panel)

    # ── the metric ───────────────────────────────────────────────────────

    def volume_density(self, orbit_coordinate: np.ndarray) -> np.ndarray:
        r"""The density :math:`c\,d\,r^{d-1}` of the chart's volume measure in the orbit coordinate.

        The derivative of the one definition
        :math:`m = c\,(T(r_{j+1}) - T(r_j))`, :math:`T(r) = r^d`.
        """
        coord = self.regions.chart.coord
        d = coord.measure_coordinate.exponent
        return coord.measure_constant * d * np.asarray(orbit_coordinate, dtype=float) ** (d - 1)

    @cached_property
    def mass(self) -> np.ndarray:
        r"""The Gram matrix :math:`W_{ij} = \int u_i u_j\,\mathrm{d}V`, ``(N, N)``, block-diagonal by panel."""
        rule = composite_gauss_legendre(self.partition.breakpoints, self.per_panel + 1)
        panel = np.repeat(np.arange(self.n_panels), self.per_panel + 1)
        shape = (self.n_panels, self.per_panel + 1)
        u = self.values(rule.pts, panel).reshape(*shape, self.per_panel)
        weight = (rule.wts * self.volume_density(rule.pts)).reshape(shape)
        return block_diag(*np.einsum("Pqi,Pq,Pqj->Pij", u, weight, u))


__all__ = ["PanelBasis"]
