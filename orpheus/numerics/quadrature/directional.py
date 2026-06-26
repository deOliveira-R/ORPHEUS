r"""Canonical directional :class:`Quadrature` for SN transport.

ONE class — :class:`Quadrature` — wrapping a :class:`DiscreteMeasure`
(the mathematical primitive) with the SN-side derived data used by
the sweep, the BC realisers, and the harmonic-moment field:

* **Reflection partners** — for each direction-cosine axis, the
  partner-ordinate index under reflection. Cached at construction
  (the nearest-neighbour search for sphere quadratures is O(N²);
  paying it once per quadrature is cheap, paying it per sweep call
  is not).
* **Octant partition** — disjoint sign-of-direction decomposition
  of the sphere cubature. Cached lazily via the underlying
  :meth:`DiscreteMeasure.partition_by`.
* **Spherical-harmonics evaluation** — the :math:`(N, L+1, 2L+1)`
  table of real spherical harmonics evaluated at the ordinates. The
  P_N Galerkin reconstruction and the harmonic-moment field consume
  this.
* **Level structure** — for cylindrical SN, the per-level
  :math:`(\mu_p, \mathrm{indices})` metadata returned alongside the
  measure by :func:`level_symmetric_sn` / :func:`product_mu_phi`.

The four families of SN quadrature (GL1D / Lebedev / Level-Symmetric
/ Product) become **classmethod factories** on a single
:class:`Quadrature` class — there is no per-family Python subclass
(R-1 Phase A detour-C: full retirement of the SN-side adapter
hierarchy). Construction sites change from
``GaussLegendre1D.create(8)`` to ``Quadrature.gauss_legendre(8)`` and
similar.

Why a single class instead of free functions
--------------------------------------------

The :class:`DiscreteMeasure` primitive is enough mathematically; the
class only adds:

1. A *cache home* for reflection partners (O(N²) to compute on
   sphere quadratures, called O(1) times per ordinate during a sweep).
2. A *bundling* for the cylindrical-level side-channel
   (:class:`LevelStructure`) that one of two quadrature families
   carries.
3. Convenience methods (``axis_cosines``, ``spherical_harmonics``)
   that take an axis index or :math:`L` and return the right slice
   of derived data.

These are all session-time concerns rather than primitive-time math.
Keeping them on a thin class is structurally honest: one cache, one
bundle, one access surface.

Legacy denormalised views
-------------------------

The :class:`Quadrature` class keeps ``mu_x`` / ``mu_y`` / ``mu_z`` /
``weights`` / ``N`` as ``@property`` views over the underlying
``measure`` and ``measure.nodes`` columns. These names are SN slab
convention; for cylindrical SN they are mis-leading (LS column 0 is
the *radial* cosine :math:`\eta`, NOT a Cartesian X-projection).
They survive only because the legacy consumer call sites (~150 reads
across orpheus/+tests/) use these names. Migration to
:meth:`axis_cosines` is a follow-up. The Pattern 7 violation that
made them dangerous (denormalised dataclass *fields* in four
parallel adapter classes) is closed by construction here: there is
exactly ONE producer of the ordinate data, the ``measure``; the
named accessors are derived views with no separate storage.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np

from orpheus.numerics.basis.spherical_harmonic_basis import (
    SphericalHarmonicBasis,
)
from orpheus.numerics.measure import (
    SPACE_SPHERE,
    DiscreteMeasure,
    DiscreteMeasurePartition,
)

if TYPE_CHECKING:
    from orpheus.numerics.frame import GalerkinFrame

from .rules_1d import gauss_legendre_on_mu
from .rules_product import product_mu_phi
from .rules_sphere import LevelStructure, lebedev_sphere, level_symmetric_sn

# Threshold below which a direction-cosine component is treated as
# zero for octant labelling. Matches the pure-z degenerate-ordinate
# threshold in ``orpheus.sn.sweep`` and ``orpheus.sn.spatial.diamond``
# (``_DEGENERATE_ABS_MU_THRESHOLD``); keep in lockstep.
_OCTANT_SIGN_EPS = 1e-15


# ═══════════════════════════════════════════════════════════════════════
# Internal helpers
# ═══════════════════════════════════════════════════════════════════════


def _octant_sign_predicate(nodes: np.ndarray) -> np.ndarray:
    r"""Sign-of-direction labelling predicate for octant partitioning.

    Returns ``+1`` for components ``> +eps``, ``-1`` for components
    ``< -eps``, and ``0`` for ``|component| <= eps``. The result has
    the same shape as ``nodes`` (``(N,)`` for 1-D slab quadratures,
    ``(N, 3)`` for spherical cubatures), so it is consumed by both
    branches of :meth:`DiscreteMeasure.partition_by`.

    Pure-axis ordinates (e.g. Lebedev's :math:`(0, 0, \pm 1)`) get a
    partition label with one or more zero components, separating them
    from the eight full octants — this is exactly the structure the
    2-D wavefront sweep relies on to short-circuit the streaming step
    for ordinates with :math:`|\mu_x| < \epsilon \wedge |\mu_y| <
    \epsilon`.
    """
    out = np.zeros(nodes.shape, dtype=int)
    out[nodes > _OCTANT_SIGN_EPS] = 1
    out[nodes < -_OCTANT_SIGN_EPS] = -1
    return out


def _find_reflections(
    tx: np.ndarray, ty: np.ndarray, tz: np.ndarray,
    rx: np.ndarray, ry: np.ndarray, rz: np.ndarray,
) -> np.ndarray:
    """Find index of closest match in ``(rx, ry, rz)`` for each
    ``(tx, ty, tz)``.

    Used by sphere-quadrature factories (Lebedev / level-symmetric /
    product) to precompute reflection partners
    (``ref[n] = arg min_k ||R x_n - x_k||`` for the requested
    reflection :math:`R`). Conceptually the same operation as a
    :meth:`DiscreteMeasure.pushforward` followed by a snapping step;
    the explicit nearest-neighbour search lives here in the
    quadrature layer because the SN consumer requires *integer*
    indices into the original node array, not a new measure with
    permuted nodes.
    """
    n = len(tx)
    ref = np.empty(n, dtype=int)
    for i in range(n):
        dist = (rx - tx[i]) ** 2 + (ry - ty[i]) ** 2 + (rz - tz[i]) ** 2
        ref[i] = np.argmin(dist)
    return ref


def _resolve_axis_to_index(axis: int | str) -> int:
    """Convert legacy ``'x'``/``'y'``/``'z'`` to ``0``/``1``/``2``.

    The :meth:`Quadrature.reflection_index` accessor accepts either
    a column index (dim-agnostic, preferred for new code) or the
    legacy SN slab tag (back-compat for the unmigrated sweep paths).
    """
    if isinstance(axis, (int, np.integer)):
        return int(axis)
    if axis == "x":
        return 0
    if axis == "y":
        return 1
    if axis == "z":
        return 2
    raise ValueError(
        f"Unknown axis label {axis!r}; expected 'x'/'y'/'z' or int 0/1/2."
    )


def _compute_sphere_reflection_partners(
    nodes: np.ndarray,
) -> dict[int, np.ndarray]:
    """Compute reflection-partner indices for axes 0/1/2 of a 3-D
    sphere cubature.

    Returns a dict keyed by axis index. Used by the sphere-quadrature
    factories (Lebedev / level-symmetric / product) to eager-populate
    the :attr:`Quadrature.reflection_partners` cache at construction.
    """
    mu_x = nodes[:, 0]
    mu_y = nodes[:, 1]
    mu_z = nodes[:, 2]
    return {
        0: _find_reflections(-mu_x, mu_y, mu_z, mu_x, mu_y, mu_z),
        1: _find_reflections(mu_x, -mu_y, mu_z, mu_x, mu_y, mu_z),
        2: _find_reflections(mu_x, mu_y, -mu_z, mu_x, mu_y, mu_z),
    }


# ═══════════════════════════════════════════════════════════════════════
# Quadrature — the single directional-quadrature class
# ═══════════════════════════════════════════════════════════════════════


@dataclass
class Quadrature:
    r"""Directional quadrature: a :class:`DiscreteMeasure` on the
    sphere or polar interval, with cached SN-derived data.

    Construction is via the named classmethod factories:

    * :meth:`gauss_legendre` — slab 1-D polar quadrature
      (Gauss-Legendre on :math:`[-1, 1]`).
    * :meth:`lebedev` — 3-D :math:`O_h`-invariant sphere cubature.
    * :meth:`level_symmetric` — Carlson-Lathrop S_N triangular
      cubature (cylindrical SN compatible).
    * :meth:`product` — Gauss-Legendre :math:`(\mu)`
      :math:`\times` equispaced :math:`(\phi)` (cylindrical SN
      compatible).

    Direct construction (``Quadrature(measure=..., ...)``) is
    permitted for cases where a custom :class:`DiscreteMeasure` is
    being wrapped (e.g. derived measures from
    :meth:`DiscreteMeasure.restrict` or
    :meth:`DiscreteMeasure.pushforward`).

    Attributes
    ----------
    measure : DiscreteMeasure
        The canonical mathematical data. Carries ``nodes``,
        ``weights``, ``space`` (``"[-1,1]"`` or ``"S^2"``), and
        polynomial-exactness metadata. **Single source of truth** —
        every other attribute on :class:`Quadrature` derives from
        this.
    reflection_partners : dict[int, np.ndarray]
        Maps axis index to the partner-ordinate index array under
        reflection across that axis. Populated by the factory at
        construction (eager precomputation pays the O(N²)
        nearest-neighbour cost once instead of per sweep call).
    level_structure : LevelStructure | None
        Per-level metadata for the cylindrical SN sweep. ``None``
        for slab / 2-D quadratures.
    """

    measure: DiscreteMeasure
    reflection_partners: dict[int, np.ndarray] = field(default_factory=dict)
    level_structure: LevelStructure | None = None

    # ────────────────────────────────────────────────────────────
    # Canonical scalar / array accessors
    # ────────────────────────────────────────────────────────────

    @property
    def n_ordinates(self) -> int:
        """Number of ordinates :math:`N`."""
        return self.measure.n_points

    @property
    def N(self) -> int:  # noqa: N802  — SN convention
        """Alias for :attr:`n_ordinates`. Matches the :math:`N` of
        :math:`S_N` nomenclature; kept for the load-bearing 1-letter
        reading at hot-loop call sites."""
        return self.measure.n_points

    @property
    def weights(self) -> np.ndarray:
        r"""Per-ordinate quadrature weights, shape ``(N,)``."""
        return self.measure.weights

    @property
    def nodes(self) -> np.ndarray:
        r"""Per-ordinate direction components, shape ``(N,)`` (1-D
        slab) or ``(N, d)`` (multi-dim sphere)."""
        return self.measure.nodes

    @property
    def dim(self) -> int:
        r"""Direction-space dimensionality :math:`d`. Returns ``1``
        for slab quadratures (1-D scalar nodes), ``3`` for sphere
        cubatures."""
        return self.measure.dim

    def axis_cosines(self, axis_index: int) -> np.ndarray:
        r"""Direction cosines along axis :math:`i` (dim-agnostic).

        For 1-D slab quadratures (``nodes.shape == (N,)``), only
        :math:`i = 0` is meaningful; other indices return zeros. For
        multi-dim cubatures, returns ``nodes[:, axis_index]`` with
        a zero fallback for axes beyond the measure's intrinsic
        dimensionality (the dim-agnostic shape primitive in
        :mod:`orpheus.transport.mesh.axis` interprets "no quadrature data on
        this axis" as "no ordinate is outflowing on it").

        This is the **canonical** per-axis accessor for new
        dim-agnostic code. The legacy :attr:`mu_x` / :attr:`mu_y` /
        :attr:`mu_z` properties are back-compat views over this
        function.
        """
        nodes = self.measure.nodes
        n_points = self.measure.n_points
        if nodes.ndim == 1:
            return nodes if axis_index == 0 else np.zeros(n_points)
        if axis_index >= nodes.shape[1]:
            return np.zeros(n_points)
        return nodes[:, axis_index]

    # ────────────────────────────────────────────────────────────
    # Legacy mu_x / mu_y / mu_z (back-compat @property views)
    # ────────────────────────────────────────────────────────────

    @property
    def mu_x(self) -> np.ndarray:
        r"""Axis-0 direction cosines. Legacy SN slab convention name.

        Equivalent to :meth:`axis_cosines(0)`. For cylindrical SN
        this is :math:`\eta` (the radial cosine), NOT a Cartesian-X
        projection — the column index, not the name, is the actual
        semantic.
        """
        return self.axis_cosines(0)

    @property
    def mu_y(self) -> np.ndarray:
        r"""Axis-1 direction cosines. Legacy SN slab convention name."""
        return self.axis_cosines(1)

    @property
    def mu_z(self) -> np.ndarray:
        r"""Axis-2 direction cosines. Legacy SN slab convention name.

        For cylindrical SN this is the axial cosine :math:`\mu` —
        the Cartesian-Z label aligns with the cylindrical convention
        here (axial ≡ z), so the name is not misleading.
        """
        return self.axis_cosines(2)

    # ────────────────────────────────────────────────────────────
    # Cylindrical-SN frame aliases (Bailey 2009 / Hébert convention)
    # ────────────────────────────────────────────────────────────
    #
    # In cylindrical SN, the natural direction-cosine names are
    # (η, ξ, μ) — radial, azimuthal, axial. These are the SAME
    # three columns of the underlying ``measure.nodes`` as
    # (mu_x, mu_y, mu_z), but the cylindrical-frame names are
    # honest about what they represent:
    #
    #   η = Ω · r̂  (radial cosine)       = nodes[:, 0]
    #   ξ = Ω · φ̂  (azimuthal cosine)    = nodes[:, 1]
    #   μ = Ω · ẑ   (axial cosine)        = nodes[:, 2]  ← same as mu_z
    #
    # The ``mu_x`` name is **misleading in cylindrical context** —
    # column 0 is the radial cosine, NOT a Cartesian-X projection.
    # Cylindrical-context call sites should read ``quad.eta`` /
    # ``quad.xi`` so the consumer code is self-documenting about
    # which coordinate frame it is operating in.

    @property
    def eta(self) -> np.ndarray:
        r"""Cylindrical-frame radial cosine :math:`\eta = \Omega \cdot \hat{r}`.

        Same data as :attr:`mu_x` (column 0 of ``measure.nodes``);
        the cylindrical-frame name is honest about the geometric
        meaning. Use this in cylindrical-SN call sites where the
        slab name ``mu_x`` would mislead a reader into thinking
        the column is a Cartesian-X projection.
        """
        return self.axis_cosines(0)

    @property
    def xi(self) -> np.ndarray:
        r"""Cylindrical-frame azimuthal cosine :math:`\xi = \Omega \cdot \hat{\varphi}`.

        Same data as :attr:`mu_y` (column 1 of ``measure.nodes``);
        the cylindrical-frame name is honest about the geometric
        meaning.
        """
        return self.axis_cosines(1)

    # ────────────────────────────────────────────────────────────
    # Reflection partners
    # ────────────────────────────────────────────────────────────

    def reflection_index(self, axis: int | str) -> np.ndarray:
        r"""Partner-ordinate index under reflection across ``axis``.

        ``axis`` accepts either a column index (``0`` / ``1`` /
        ``2`` — dim-agnostic, preferred for new code) or the legacy
        SN slab tag (``"x"`` / ``"y"`` / ``"z"`` — back-compat for
        the unmigrated sweep paths).

        For sphere cubatures (Lebedev / level-symmetric / product)
        the partner indices are precomputed at construction via
        nearest-neighbour matching against the reflected node
        positions. For slab GL1D the partners are derived from the
        :math:`i \leftrightarrow N - 1 - i` symmetry of the GL nodes.
        """
        axis_idx = _resolve_axis_to_index(axis)
        try:
            return self.reflection_partners[axis_idx]
        except KeyError as exc:
            raise ValueError(
                f"Quadrature has no precomputed reflection partner for "
                f"axis index {axis_idx}; available axes: "
                f"{sorted(self.reflection_partners.keys())}."
            ) from exc

    # ────────────────────────────────────────────────────────────
    # Spherical harmonics evaluation
    # ────────────────────────────────────────────────────────────

    def spherical_harmonics(self, L: int) -> np.ndarray:
        r"""Real spherical harmonics :math:`Y_l^m` evaluated at the
        ordinates, shape :math:`(N, L+1, 2L+1)`.

        For slab GL1D quadratures only the :math:`m = 0` harmonics
        :math:`Y_l^0(\mu_x)` carry non-zero values; the other slots
        are filled with zeros (1-D angular variation is purely
        polar). For sphere cubatures all :math:`(l, m)` pairs are
        evaluated.
        """
        return SphericalHarmonicBasis(L=L).evaluate_from_components(
            self.axis_cosines(0),
            self.axis_cosines(1),
            self.axis_cosines(2),
        )

    def angular_frame(self, L: int) -> "GalerkinFrame":
        r"""The degree-:math:`L` spherical-harmonic :class:`~orpheus.numerics.frame.GalerkinFrame` on this quadrature.

        Binds :class:`~orpheus.numerics.basis.SphericalHarmonicBasis` (degree
        :math:`L`) to this quadrature's :math:`S^2` measure: the ``(N, 3)``
        direction cosines (a slab's polar :math:`\mu` embeds as
        :math:`(\mu, 0, 0)` — the SAME column-stacked embedding
        :meth:`spherical_harmonics` uses internally) carrying the quadrature
        weights as the analysis metric. The frame is a pure-Galerkin frame
        (test IS trial — the SH basis is its own test basis); its
        :attr:`~orpheus.numerics.frame.FrameBase.analysis` face is the
        :math:`Y^* W` moment projection :math:`M` and its
        :attr:`~orpheus.numerics.frame.FrameBase.reconstruction` face the
        addition-theorem synthesis :math:`R`; ``frame.table`` equals
        :meth:`spherical_harmonics` ``(L)`` bit-identically (both route
        through :meth:`SphericalHarmonicBasis.evaluate` on these cosines).

        The single source of the angular frame consumed by
        :class:`~orpheus.transport.operators.scattering.ScatteringOperator` — its §5.6 kernel
        AND the in-sweep moment accumulation share THIS object, so the
        projection table is never re-evaluated or allowed to diverge.

        Naming tripwire: ``angular_frame`` names the quadrature's permanent
        physical **axis** (the :math:`S^2` directions), NOT the contingent
        spherical-harmonic basis — preserving the greppable
        ``angular_frame``/``spatial_frame``/``energy_frame`` family. The
        :class:`SphericalHarmonicBasis` hardcode is honest while the angular
        axis has exactly one basis; a *second* angular basis (a Fourier /
        Slepian angular expansion) parametrises the basis —
        ``angular_frame(basis=...)``, a SIGNATURE change — and must NOT rename
        the method.
        """
        from orpheus.numerics.frame import GalerkinFrame

        s2_measure = DiscreteMeasure(
            nodes=np.column_stack(
                [self.axis_cosines(0), self.axis_cosines(1), self.axis_cosines(2)]
            ),
            weights=self.weights,
            support=SPACE_SPHERE,
        )
        return GalerkinFrame(SphericalHarmonicBasis(L=L), s2_measure)

    # ────────────────────────────────────────────────────────────
    # Octants (cached partition by sign-of-direction)
    # ────────────────────────────────────────────────────────────

    @cached_property
    def octants(self) -> tuple[DiscreteMeasurePartition, ...]:
        r"""Disjoint partition of this quadrature by sign-of-direction.

        Realises the direct-sum decomposition

        .. math::

            \mu_{S^2} \;=\; \bigoplus_{\sigma} \mu_\sigma,
            \qquad \sigma = (\mathrm{sign}\,\mu_x,
                              \mathrm{sign}\,\mu_y,
                              \mathrm{sign}\,\mu_z)

        Pure-axis ordinates with :math:`|\mu_{axis}| < 10^{-15}`
        carry a ``0`` component in their octant label and form their
        own partition entry. Cached lazily — partitioning is
        mesh-time work, paid once.
        """
        return self.measure.partition_by(_octant_sign_predicate)

    # ────────────────────────────────────────────────────────────
    # Level-structure passthroughs (cylindrical SN)
    # ────────────────────────────────────────────────────────────

    @property
    def n_levels(self) -> int:
        r"""Number of polar levels (cylindrical SN side-channel).

        ``1`` for slab quadratures (the GL :math:`\mu` axis is the
        only "level"). For cylindrical-compatible cubatures
        (level-symmetric, product) returns
        :attr:`LevelStructure.n_levels`.
        """
        if self.level_structure is None:
            return 1
        return self.level_structure.n_levels

    @property
    def level_indices(self) -> list[np.ndarray]:
        r"""Per-level ordinate-index lists (cylindrical SN side-channel).

        ``[arange(N)]`` for slab quadratures (single "level"). For
        cylindrical-compatible cubatures returns
        :attr:`LevelStructure.level_indices`.
        """
        if self.level_structure is None:
            return [np.arange(self.measure.n_points)]
        return self.level_structure.level_indices

    @property
    def level_mu(self) -> np.ndarray:
        r"""Per-level polar cosines (cylindrical SN side-channel).

        ``np.array([0.0])`` for slab quadratures. For
        cylindrical-compatible cubatures returns
        :attr:`LevelStructure.level_mu`.
        """
        if self.level_structure is None:
            return np.array([0.0])
        return self.level_structure.level_mu

    # ────────────────────────────────────────────────────────────
    # Classmethod factories
    # ────────────────────────────────────────────────────────────

    @classmethod
    def gauss_legendre(cls, n_ordinates: int = 16) -> "Quadrature":
        r""":math:`N`-point Gauss-Legendre quadrature on :math:`[-1, 1]`.

        Slab transport polar quadrature. ``N`` must be even for SN
        (half-range integration over each hemisphere).

        The partner-under-x-reflection of ordinate :math:`i` is
        :math:`N - 1 - i` by GL-node symmetry; the y- and z-axes
        carry zero cosines (1-D), so every ordinate is its own
        partner under those reflections.
        """
        measure = gauss_legendre_on_mu(n_ordinates)
        N = measure.n_points
        identity = np.arange(N)
        partners = {
            0: identity[::-1].copy(),  # GL x-reflection: i ↔ N-1-i
            1: identity,               # 1-D: mu_y == 0
            2: identity,               # 1-D: mu_z == 0
        }
        return cls(
            measure=measure,
            reflection_partners=partners,
            level_structure=None,
        )

    @classmethod
    def lebedev(cls, order: int = 17) -> "Quadrature":
        r"""Lebedev :math:`O_h`-invariant cubature on the sphere
        :math:`S^2`.

        ``order`` is the polynomial-exactness order. Common values
        are 5, 7, 11, 17, 21, ... — see
        :func:`orpheus.numerics.quadrature.lebedev_sphere` for the
        full list of supported orders.

        Weights sum to :math:`4\pi`.
        """
        measure = lebedev_sphere(order)
        partners = _compute_sphere_reflection_partners(measure.nodes)
        return cls(
            measure=measure,
            reflection_partners=partners,
            level_structure=None,
        )

    @classmethod
    def level_symmetric(cls, sn_order: int = 4) -> "Quadrature":
        r"""Carlson-Lathrop level-symmetric :math:`S_N` cubature.

        ``sn_order`` (even :math:`\ge 2`) is the standard SN order;
        for :math:`S_N` the cubature has :math:`N (N + 2)`
        ordinates. Common values: 4, 8, 16.

        Provides the per-level :math:`(\mu_p, \mathrm{indices})`
        side-channel required by the cylindrical SN sweep.
        :math:`O_h`-invariant; weights sum to :math:`4\pi`.
        """
        measure, structure = level_symmetric_sn(sn_order)
        partners = _compute_sphere_reflection_partners(measure.nodes)
        return cls(
            measure=measure,
            reflection_partners=partners,
            level_structure=structure,
        )

    @classmethod
    def product(cls, n_mu: int = 8, n_phi: int = 8) -> "Quadrature":
        r"""Product quadrature: Gauss-Legendre on :math:`\mu` and
        equispaced on :math:`\phi`.

        The polar angle :math:`\theta` is discretised via GL on
        :math:`\mu = \cos\theta \in [-1, 1]`; the azimuthal angle
        :math:`\phi` is discretised uniformly on :math:`[0, 2\pi)`.

        Direction cosines (cylindrical SN convention):

        * :math:`\eta = \mu_x = \sin\theta \cos\phi` — radial cosine
        * :math:`\xi  = \mu_y = \sin\theta \sin\phi` — azimuthal cosine
        * :math:`\mu  = \mu_z = \cos\theta`           — axial cosine

        Weights sum to :math:`4\pi`. Provides the per-level
        side-channel required by the cylindrical SN sweep.
        """
        measure, structure = product_mu_phi(n_mu, n_phi)
        partners = _compute_sphere_reflection_partners(measure.nodes)
        return cls(
            measure=measure,
            reflection_partners=partners,
            level_structure=structure,
        )


__all__ = ["Quadrature"]
