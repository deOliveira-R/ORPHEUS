r"""Angular quadrature for SN transport — adapters over
:class:`~orpheus.numerics.measure.DiscreteMeasure`.

Provides the :class:`AngularQuadrature` Protocol plus four concrete
adapters:

* :class:`GaussLegendre1D` — for 1-D slab problems
  (:math:`\mu \in [-1, 1]`, weights sum to 2).
* :class:`LebedevSphere` — for 2-D / 3-D problems on
  :math:`O_h`-invariant Lebedev grids (weights sum to :math:`4\pi`).
* :class:`LevelSymmetricSN` — Carlson-Lathrop level-symmetric
  :math:`S_N` on :math:`S^2` (weights sum to :math:`4\pi`).
* :class:`ProductQuadrature` — product rule
  :math:`\mu_{\text{GL}} \times \phi_{\text{equispaced}}` on
  :math:`S^2`.

Why adapters?
-------------

The four classes pre-date the
:class:`~orpheus.numerics.measure.DiscreteMeasure` abstraction. SN
solvers, sweeps, BiCGSTAB operators, and meshes consume them through
~50 attribute-access sites in :mod:`orpheus.sn` (``quad.mu_x``,
``quad.weights``, ``quad.reflection_index('x')``, …). Issue 4 of
``.claude/plans/sn_reshape.md`` re-platforms each class as a *thin
adapter* over a measure-returning rule function from
:mod:`orpheus.numerics.quadrature`, caching numpy views on construction
so the legacy attribute API stays bit-identical (verified by the 11
regression snapshots at ``tests/sn/regression/snapshots/``).

The adapter caches in ``__init__``:

* ``mu_x`` / ``mu_y`` / ``mu_z`` — column views into
  ``measure.nodes`` (no copy).
* ``weights`` / ``N`` — view of ``measure.weights`` and
  ``measure.n_points``.
* ``_ref_x`` / ``_ref_y`` / ``_ref_z`` — reflection-index arrays via
  :meth:`DiscreteMeasure.pushforward` plus nearest-neighbour matching.
* ``level_indices`` / ``level_mu`` / ``n_levels`` — the
  :class:`~orpheus.numerics.quadrature.rules_sphere.LevelStructure`
  side-channel, for the cylindrical SN sweep.

The underlying :class:`DiscreteMeasure` is exposed as ``self.measure``
for callers that want the structurally-richer object (composability,
``invariance_group``, ``degree_of_exactness``).

The solver uses ``1/sum(weights)`` as the isotropic normalisation
factor, so it remains quadrature-agnostic.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Protocol, runtime_checkable

import numpy as np

from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.quadrature import (
    gauss_legendre_on_mu,
    lebedev_sphere,
    level_symmetric_sn,
    product_mu_phi,
)
from orpheus.numerics.quadrature.rules_sphere import LevelStructure
from orpheus.numerics.spherical_harmonics import evaluate_real_sh as _build_spherical_harmonics


@runtime_checkable
class AngularQuadrature(Protocol):
    """Contract for any angular quadrature usable by the SN solver."""

    mu_x: np.ndarray       # (N,) x-direction cosines
    mu_y: np.ndarray       # (N,) y-direction cosines (0 for 1D)
    weights: np.ndarray    # (N,) quadrature weights
    N: int                 # number of ordinates

    def reflection_index(self, axis: str) -> np.ndarray:
        """Index array: ref[n] = partner of ordinate n reflected in ``axis``."""
        ...

    def spherical_harmonics(self, L: int) -> np.ndarray:
        """(N, L+1, 2L+1) real spherical harmonics Y[n, l, l+m]."""
        ...


# ═══════════════════════════════════════════════════════════════════════
# Utility: nearest-neighbour reflection-index lookup
# ═══════════════════════════════════════════════════════════════════════


def _find_reflections(
    tx: np.ndarray, ty: np.ndarray, tz: np.ndarray,
    rx: np.ndarray, ry: np.ndarray, rz: np.ndarray,
) -> np.ndarray:
    """Find index of closest match in (rx,ry,rz) for each (tx,ty,tz).

    Used by Lebedev / level-symmetric / product adapters to precompute
    reflection partners (``ref[n] = arg min_k |R x_n - x_k|`` for the
    requested reflection ``R``). Conceptually the same operation as a
    :meth:`DiscreteMeasure.pushforward` followed by a snapping step;
    we keep the explicit nearest-neighbour search in the adapter
    layer because the SN consumer requires *integer* indices into
    the original node array, not a new measure with permuted nodes.
    """
    n = len(tx)
    ref = np.empty(n, dtype=int)
    for i in range(n):
        dist = (rx - tx[i])**2 + (ry - ty[i])**2 + (rz - tz[i])**2
        ref[i] = np.argmin(dist)
    return ref


# ═══════════════════════════════════════════════════════════════════════
# Gauss-Legendre 1-D adapter
# ═══════════════════════════════════════════════════════════════════════


@dataclass
class GaussLegendre1D:
    """Gauss-Legendre quadrature on [-1, 1] for 1D slab transport.

    Adapter over
    :func:`orpheus.numerics.quadrature.gauss_legendre_on_mu`. Caches
    ``mu_x`` (the GL nodes), ``mu_y`` (zeros), and ``weights`` as
    array views of the underlying :class:`DiscreteMeasure`. ``N`` is
    the node count. Weights sum to 2.
    """

    mu_x: np.ndarray
    mu_y: np.ndarray
    weights: np.ndarray
    N: int
    measure: DiscreteMeasure | None = field(default=None, repr=False)

    @property
    def mu(self) -> np.ndarray:
        """Alias for mu_x (the 1D direction cosines)."""
        return self.mu_x

    @classmethod
    def create(cls, n_ordinates: int = 16) -> GaussLegendre1D:
        """Build N-point GL quadrature (must be even for SN)."""
        measure = gauss_legendre_on_mu(n_ordinates)
        # ``measure.nodes`` is shape (N,) — the polar cosine. We cache
        # that as ``mu_x`` (slab convention: x is the streaming axis).
        # ``mu_y`` is identically zero in 1-D.
        return cls(
            mu_x=measure.nodes,
            mu_y=np.zeros(measure.n_points),
            weights=measure.weights,
            N=measure.n_points,
            measure=measure,
        )

    def reflection_index(self, axis: str) -> np.ndarray:
        """GL is symmetric: partner of i is N-1-i for x-reflection."""
        if axis == "x":
            return np.arange(self.N)[::-1].copy()
        else:
            # y-reflection: mu_y=0, so every ordinate is its own partner
            return np.arange(self.N)

    def spherical_harmonics(self, L: int) -> np.ndarray:
        """1D harmonics: Y_0^0=1, Y_1^0=μ_x (only x-component in 1D)."""
        return _build_spherical_harmonics(
            L, self.mu_x, self.mu_y, np.zeros(self.N),
        )


# ═══════════════════════════════════════════════════════════════════════
# Lebedev sphere adapter
# ═══════════════════════════════════════════════════════════════════════


@dataclass
class LebedevSphere:
    """Lebedev quadrature on the unit sphere for 2D/3D transport.

    Adapter over :func:`orpheus.numerics.quadrature.lebedev_sphere`.
    Caches ``mu_x`` / ``mu_y`` / ``mu_z`` as column views of the
    underlying :class:`DiscreteMeasure`'s ``nodes`` array (which is
    shape ``(N, 3)``). Weights sum to :math:`4\\pi`.
    """

    mu_x: np.ndarray
    mu_y: np.ndarray
    mu_z: np.ndarray
    weights: np.ndarray
    N: int
    _ref_x: np.ndarray
    _ref_y: np.ndarray
    _ref_z: np.ndarray
    measure: DiscreteMeasure | None = field(default=None, repr=False)

    @classmethod
    def create(cls, order: int = 17) -> LebedevSphere:
        """Build Lebedev quadrature of given polynomial-exact order."""
        measure = lebedev_sphere(order)
        nodes = measure.nodes  # (N, 3)
        mu_x = nodes[:, 0]
        mu_y = nodes[:, 1]
        mu_z = nodes[:, 2]
        n_pts = measure.n_points

        ref_x = _find_reflections(-mu_x, mu_y, mu_z, mu_x, mu_y, mu_z)
        ref_y = _find_reflections(mu_x, -mu_y, mu_z, mu_x, mu_y, mu_z)
        ref_z = _find_reflections(mu_x, mu_y, -mu_z, mu_x, mu_y, mu_z)

        return cls(
            mu_x=mu_x, mu_y=mu_y, mu_z=mu_z,
            weights=measure.weights, N=n_pts,
            _ref_x=ref_x, _ref_y=ref_y, _ref_z=ref_z,
            measure=measure,
        )

    def reflection_index(self, axis: str) -> np.ndarray:
        if axis == "x":
            return self._ref_x
        elif axis == "y":
            return self._ref_y
        elif axis == "z":
            return self._ref_z
        raise ValueError(f"Unknown axis: {axis}")

    def spherical_harmonics(self, L: int) -> np.ndarray:
        """(N, L+1, 2L+1) real spherical harmonics for Lebedev ordinates."""
        return _build_spherical_harmonics(L, self.mu_x, self.mu_y, self.mu_z)


# ═══════════════════════════════════════════════════════════════════════
# Level-Symmetric S_N adapter
# ═══════════════════════════════════════════════════════════════════════


@dataclass
class LevelSymmetricSN:
    """Level-symmetric :math:`S_N` quadrature on the unit sphere.

    Adapter over
    :func:`orpheus.numerics.quadrature.level_symmetric_sn`.
    Standard triangular quadrature with :math:`N/2` :math:`\\mu`-levels
    per hemisphere; provides the ``level_indices`` structure required
    by the cylindrical SN sweep for azimuthal redistribution
    (Bailey et al. 2009, Eq. 50). Weights sum to :math:`4\\pi`.
    """

    mu_x: np.ndarray       # η — radial direction cosines
    mu_y: np.ndarray       # ξ — azimuthal direction cosines
    mu_z: np.ndarray       # μ — axial direction cosines
    weights: np.ndarray
    N: int
    _ref_x: np.ndarray
    _ref_y: np.ndarray
    _ref_z: np.ndarray

    # Level structure
    n_levels: int
    level_indices: list[np.ndarray]
    level_mu: np.ndarray
    measure: DiscreteMeasure | None = field(default=None, repr=False)

    @classmethod
    def create(cls, sn_order: int = 4) -> LevelSymmetricSN:
        """Build :math:`S_N` level-symmetric quadrature of given order."""
        measure, structure = level_symmetric_sn(sn_order)
        nodes = measure.nodes  # (N, 3)
        mu_x = nodes[:, 0]
        mu_y = nodes[:, 1]
        mu_z = nodes[:, 2]
        N = measure.n_points

        ref_x = _find_reflections(-mu_x, mu_y, mu_z, mu_x, mu_y, mu_z)
        ref_y = _find_reflections(mu_x, -mu_y, mu_z, mu_x, mu_y, mu_z)
        ref_z = _find_reflections(mu_x, mu_y, -mu_z, mu_x, mu_y, mu_z)

        return cls(
            mu_x=mu_x, mu_y=mu_y, mu_z=mu_z,
            weights=measure.weights, N=N,
            _ref_x=ref_x, _ref_y=ref_y, _ref_z=ref_z,
            n_levels=structure.n_levels,
            level_indices=structure.level_indices,
            level_mu=structure.level_mu,
            measure=measure,
        )

    def reflection_index(self, axis: str) -> np.ndarray:
        if axis == "x":
            return self._ref_x
        elif axis == "y":
            return self._ref_y
        elif axis == "z":
            return self._ref_z
        raise ValueError(f"Unknown axis: {axis}")

    def spherical_harmonics(self, L: int) -> np.ndarray:
        return _build_spherical_harmonics(L, self.mu_x, self.mu_y, self.mu_z)


# ═══════════════════════════════════════════════════════════════════════
# Product Quadrature (GL in μ × equispaced in φ) adapter
# ═══════════════════════════════════════════════════════════════════════


@dataclass
class ProductQuadrature:
    r"""Product quadrature: Gauss-Legendre :math:`(\mu)` :math:`\times`
    equispaced :math:`(\phi)`.

    Adapter over :func:`orpheus.numerics.quadrature.product_mu_phi`.
    The polar angle :math:`\theta` is discretised via Gauss-Legendre
    on :math:`\mu = \cos\theta \in [-1, 1]`; the azimuthal angle
    :math:`\phi` is discretised uniformly on :math:`[0, 2\pi)`.

    Direction cosines:

    * :math:`\mu_z = \mu` (axial, :math:`= \cos\theta`)
    * :math:`\mu_x = \eta = \sin\theta\cos\phi` (radial for cylindrical)
    * :math:`\mu_y = \xi = \sin\theta\sin\phi` (azimuthal for cylindrical)

    Weights: :math:`w = w_{\text{GL}}(\mu) \cdot (2\pi / n_\phi)` —
    sum to :math:`4\pi`.

    Provides ``level_indices`` for the cylindrical sweep.
    """

    mu_x: np.ndarray       # η = sin(θ)cos(φ)
    mu_y: np.ndarray       # ξ = sin(θ)sin(φ)
    mu_z: np.ndarray       # μ = cos(θ)
    weights: np.ndarray
    N: int
    _ref_x: np.ndarray
    _ref_y: np.ndarray
    _ref_z: np.ndarray

    # Level structure
    n_levels: int
    level_indices: list[np.ndarray]
    level_mu: np.ndarray
    measure: DiscreteMeasure | None = field(default=None, repr=False)

    @classmethod
    def create(cls, n_mu: int = 8, n_phi: int = 8) -> ProductQuadrature:
        """Build product quadrature with ``n_mu`` GL points and
        ``n_phi`` azimuthal points.

        Parameters
        ----------
        n_mu : int
            Number of Gauss-Legendre points in :math:`\\mu` (polar).
        n_phi : int
            Number of equispaced points in :math:`\\phi` (azimuthal).
        """
        measure, structure = product_mu_phi(n_mu, n_phi)
        nodes = measure.nodes  # (N_total, 3)
        mu_x = nodes[:, 0]
        mu_y = nodes[:, 1]
        mu_z = nodes[:, 2]
        n_total = measure.n_points

        ref_x = _find_reflections(-mu_x, mu_y, mu_z, mu_x, mu_y, mu_z)
        ref_y = _find_reflections(mu_x, -mu_y, mu_z, mu_x, mu_y, mu_z)
        ref_z = _find_reflections(mu_x, mu_y, -mu_z, mu_x, mu_y, mu_z)

        return cls(
            mu_x=mu_x, mu_y=mu_y, mu_z=mu_z,
            weights=measure.weights, N=n_total,
            _ref_x=ref_x, _ref_y=ref_y, _ref_z=ref_z,
            n_levels=structure.n_levels,
            level_indices=structure.level_indices,
            level_mu=structure.level_mu,
            measure=measure,
        )

    def reflection_index(self, axis: str) -> np.ndarray:
        if axis == "x":
            return self._ref_x
        elif axis == "y":
            return self._ref_y
        elif axis == "z":
            return self._ref_z
        raise ValueError(f"Unknown axis: {axis}")

    def spherical_harmonics(self, L: int) -> np.ndarray:
        return _build_spherical_harmonics(L, self.mu_x, self.mu_y, self.mu_z)
