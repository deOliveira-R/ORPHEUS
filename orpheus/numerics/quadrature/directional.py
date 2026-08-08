r"""Canonical directional :class:`Quadrature` for SN transport.

ONE class — :class:`Quadrature` — wrapping a :class:`DiscreteMeasure`
(the mathematical primitive) with the SN-side derived data used by
the sweep, the BC realisers, and the harmonic-moment field:

* **Ordinate permutations** — :meth:`Quadrature.ordinate_permutation`
  answers "which permutation does a rigid motion induce on this rule's
  weighted ordinates?", certified (every image matches a node, the
  match is a bijection, the matched weights are equal) or ``None``.
  The boundary realizer, the specular certifications, and the
  curvilinear pole seed all read this one source; consumers that pay
  it per iteration hoist the derived ``.indices`` to construction
  (e.g. ``_OneDimScanWalk._ensure_pole_mirror``). Until G6.3 §7d a
  precomputed per-axis mirror table (``reflection_partners`` /
  ``reflection_index``) carried the same answers for the three
  axis-mirrors only; the general method retired it.
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

1. The *SN-side derived-data surface* (``ordinate_permutation``,
   ``spherical_harmonics``, octants) over the one underlying measure.
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

from dataclasses import dataclass
from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np

from orpheus.numerics.basis.spherical_harmonic_basis import (
    MirrorEvenSphericalHarmonicBasis,
    SphericalHarmonicBasis,
)
from orpheus.numerics.measure import (
    SPACE_SPHERE,
    DiscreteMeasure,
    DiscreteMeasurePartition,
)
# An ordinate permutation is a group action on the measure, so it is
# certified by the SAME machinery that proves invariance — one source of
# truth for "does this motion permute these weighted nodes?"
# (``RigidMotion.preserves``), asked with the measure-level windows and
# the canonical R³ embedding (``_embedded_nodes``).
from orpheus.numerics.symmetry import (
    _NODE_WINDOW_FACTOR,
    SubgroupOfO3,
    _embedded_nodes,
)

if TYPE_CHECKING:
    from orpheus.geometry.transformation import Permutation, RigidMotion
    from orpheus.numerics.frame import GalerkinFrame

from .rules_1d import gauss_legendre_on_mu
from .rules_circle import STAGGERED, periodic_trapezoid
from .rules_product import product_mu_phi, spherical_product
from .rules_sphere import LevelStructure, lebedev_sphere, level_symmetric_sn

# Threshold below which a direction-cosine component is treated as
# zero for octant labelling. Matches the pure-z degenerate-ordinate
# threshold in ``orpheus.sn.sweep`` and ``orpheus.transport.spatial.diamond``
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


#: Node-match tolerance for certifying an ordinate permutation. Passed to
#: :meth:`~orpheus.geometry.transformation.RigidMotion.preserves`, which
#: matches positions inside ``atol`` (fed ``_REFLECTION_ATOL *
#: _NODE_WINDOW_FACTOR`` below) and compares weights at ``weight_atol``
#: (fed ``_REFLECTION_ATOL``) — the same two windows the measure-level
#: orbit certification uses.
_REFLECTION_ATOL = 1e-13


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
    level_structure : LevelStructure | None
        Per-level metadata for the cylindrical SN sweep. ``None``
        for slab / 2-D quadratures.
    """

    measure: DiscreteMeasure
    level_structure: LevelStructure | None = None
    #: The symmetry subgroup this rule was QUOTIENTED by, or ``None`` for
    #: an unfolded rule.  Set by :meth:`quotient`; consumed by the
    #: harmonic machinery (:meth:`spherical_harmonics` /
    #: :meth:`angular_frame`), which must bind the σ-EVEN sub-basis on a
    #: folded rule — the σ-odd harmonics are not in the quotient's
    #: function space and their raw moments are garbage, not zero (Q5.6).
    folded_by: "SubgroupOfO3 | None" = None

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
    # Ordinate permutations
    # ────────────────────────────────────────────────────────────

    def ordinate_permutation(
        self, motion: "RigidMotion"
    ) -> "Permutation | None":
        r"""The permutation ``motion`` induces on this rule's weighted
        ordinates, or ``None`` if the motion does not carry them onto
        themselves.

        **The single source for "which ordinate permutation does a rigid
        motion induce?"** (G6.3 step 7). Before this method the tree had two
        paths to that answer — a certified axis-mirror table
        (``reflection_index``, precomputed at construction and retired at
        §7d.3) and nothing at all for any other motion — so a translation
        deck element realized through a hand-argued identity and a sector
        rotation had no path. This method is the general question, asked
        with exactly the discipline the table carried:

        * the same canonical :math:`\mathbb{R}^3` embedding
          (:func:`~orpheus.numerics.symmetry._embedded_nodes` — a polar
          :math:`\mu` becomes :math:`(\mu, 0, 0)`), with a lower-dimensional
          motion lifted through
          :meth:`~orpheus.geometry.transformation.RigidMotion.embedded_in`;
        * the same two tolerance windows (node window
          :data:`~orpheus.numerics.symmetry._NODE_WINDOW_FACTOR` times the
          weight window, the measure-level claim recorded on
          ``_orbit_closure``);
        * the same three certifications, delivered by
          :meth:`~orpheus.geometry.transformation.RigidMotion.preserves`:
          every image matches a node (no bare ``argmin`` — ERR-074), the
          match is a **bijection** (ERR-073), and the matched nodes carry
          equal weights (a deck transformation is measure-preserving —
          ERR-042's third check).

        Returning the :class:`~orpheus.geometry.transformation.Permutation`
        rather than a ``bool`` or a bare index array is the L-013 lesson:
        the returned object PROVES its own bijectivity at construction, and
        every downstream reading (the boundary trace arrow, orbit structure,
        fixed points) is a reading of :math:`\pi`. ``None`` is the honest
        "this motion is not a symmetry of this rule" — callers own the
        domain-specific refusal message.

        Convention (pinned by the core's homomorphism law):
        :math:`g(\Omega_i) = \Omega_{\pi(i)}` — ``pi.indices[i]`` is the
        index of ordinate :math:`i`'s IMAGE.

        The match runs through
        :attr:`~orpheus.geometry.transformation.RigidMotion.linear_part`,
        because ordinates are DIRECTIONS: the translation does not act on
        them (a wrap's ordinate permutation is the identity — `[M]` feeding
        the affine element instead translates the nodes off the sphere and
        answers ``None``, a false "not a symmetry").
        """
        motion3 = (
            motion if motion.dimension == 3 else motion.embedded_in(3)
        )
        return motion3.linear_part.preserves(
            _embedded_nodes(self.measure),
            self.measure.weights,
            atol=_REFLECTION_ATOL * _NODE_WINDOW_FACTOR,
            weight_atol=_REFLECTION_ATOL,
        )

    # ────────────────────────────────────────────────────────────
    # The fold — quotient by a symmetry subgroup (Q5.6)
    # ────────────────────────────────────────────────────────────

    def quotient(
        self, group: "SubgroupOfO3", *, atol: float = 1e-13
    ) -> "Quadrature":
        r"""The quotient quadrature :math:`\mu / G` — **THE FOLD**, lifted
        to the :class:`Quadrature` tier.

        Delegates both halves to the landed primitives and adds nothing:

        * the **measure** folds by orbit-stabilizer
          (:meth:`~orpheus.numerics.measure.DiscreteMeasure.quotient` —
          total mass preserved exactly, weights *derived* as
          :math:`W = w \cdot |G| / |\mathrm{Stab}|`, the exactness claim
          dropped, the symmetry consumed);
        * the **level structure** DESCENDS by selection
          (:meth:`~orpheus.numerics.quadrature.rules_sphere.LevelStructure.quotient`
          — charts bit-copied, each level's order the parent's η-order
          RESTRICTED to the surviving representatives, fiberwise
          certified by per-level mass).

        Every refusal arm lives on those primitives (a non-invariant
        measure, a continuous group, a 1-D measure without the 3-D
        realization, a level-merging fold, a foreign measure); this
        lift owns none.  A rule without a level structure (Lebedev)
        folds its measure alone and keeps ``level_structure=None``.

        For the cylindrical σ_y-fold with the DERIVED staggered offset,
        use :meth:`folded_product` — it co-selects the offset the free
        fold requires (T25) so the incoherent NODE_ALIGNED-plus-fold
        combination is not spellable.
        """
        folded_measure = self.measure.quotient(group, atol=atol)
        folded_structure = (
            self.level_structure.quotient(
                parent=self.measure, onto=folded_measure
            )
            if self.level_structure is not None
            else None
        )
        return Quadrature(
            measure=folded_measure,
            level_structure=folded_structure,
            folded_by=group,
        )

    # ────────────────────────────────────────────────────────────
    # Spherical harmonics evaluation
    # ────────────────────────────────────────────────────────────

    def _harmonic_basis(self, L: int) -> SphericalHarmonicBasis:
        r"""The harmonic basis this rule's function space admits.

        An unfolded rule spans the full degree-:math:`L` real SH basis.
        A FOLDED rule (:attr:`folded_by` set) is a quotient measure —
        the σ-odd harmonics are not in its function space, and the raw
        :math:`Y^{\mathsf T} W` analysis (which has no Gram division
        anywhere on the scattering path) would return garbage for them,
        not zero — so the fold binds the σ-EVEN sub-basis
        (:class:`~orpheus.numerics.basis.spherical_harmonic_basis.MirrorEvenSphericalHarmonicBasis`,
        rectangular layout, odd columns structurally zeroed).  Only a
        single-mirror fold has the ±parity split this realizes; folding
        by any other group refuses here until a consumer exists.
        """
        if self.folded_by is None:
            return SphericalHarmonicBasis(L=L)
        axis = self.folded_by.mirror_axis
        if axis is None:
            raise NotImplementedError(
                f"harmonic machinery on a rule folded by "
                f"{self.folded_by.name}: only a single-mirror fold has "
                f"the even/odd parity split the restricted sub-basis "
                f"realizes; no consumer needs another group's quotient "
                f"basis yet."
            )
        return MirrorEvenSphericalHarmonicBasis(L=L, mirror_axis=axis)

    def spherical_harmonics(self, L: int) -> np.ndarray:
        r"""Real spherical harmonics :math:`Y_l^m` evaluated at the
        ordinates, shape :math:`(N, L+1, 2L+1)`.

        For slab GL1D quadratures only the :math:`m = 0` harmonics
        :math:`Y_l^0(\mu_x)` carry non-zero values; the other slots
        are filled with zeros (1-D angular variation is purely
        polar). For sphere cubatures all :math:`(l, m)` pairs are
        evaluated.
        """
        return self._harmonic_basis(L).evaluate_from_components(
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
        return GalerkinFrame(self._harmonic_basis(L), s2_measure)

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

        The x-mirror pairs ordinate :math:`i` with :math:`N - 1 - i`
        by GL-node symmetry; the y- and z-mirrors fix every ordinate
        (1-D nodes embed as :math:`(\mu, 0, 0)`). Both facts are
        DERIVED by :meth:`ordinate_permutation`, not stored.
        """
        return cls(
            measure=gauss_legendre_on_mu(n_ordinates),
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
        return cls(
            measure=lebedev_sphere(order),
            level_structure=None,
        )

    @classmethod
    def level_symmetric(cls, sn_order: int = 4) -> "Quadrature":
        r"""Carlson-Lathrop level-symmetric :math:`S_N` cubature.

        ``sn_order`` (even :math:`\ge 2`) is the standard SN order;
        for :math:`S_N` the cubature has :math:`N (N + 2)`
        ordinates. Supported: :math:`S_2` through :math:`S_{18}` —
        the family **refuses above** :math:`S_{18}` (the per-orbit
        weight solve goes negative on the moment-matched nodes from
        :math:`S_{20}`, and :math:`S_{22}` is intrinsically dead for
        the whole even-moment family; #327/#337, and see
        :func:`~orpheus.numerics.quadrature.rules_sphere.level_symmetric_sn`).

        Provides the per-level :math:`(\mu_p, \mathrm{indices})`
        side-channel required by the cylindrical SN sweep.
        :math:`O_h`-invariant; weights sum to :math:`4\pi`.
        """
        measure, structure = level_symmetric_sn(sn_order)
        return cls(measure=measure, level_structure=structure)

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
        return cls(measure=measure, level_structure=structure)

    @classmethod
    def folded_product(cls, n_mu: int = 8, n_phi: int = 8) -> "Quadrature":
        r"""The σ_y-QUOTIENT of the staggered product rule — the
        cylindrical fold.

        Builds :math:`\mathrm{GL}(n_\mu) \times
        \mathrm{trapezoid}(n_\varphi,\ \mathrm{STAGGERED})` on
        :math:`S^2` and quotients it by the azimuthal mirror
        :math:`\sigma_y` (:math:`\xi \to -\xi`).  ``n_phi`` counts the
        PARENT circle's nodes; the fold keeps one representative per
        mirror orbit, so the result carries
        :math:`n_\mu \cdot n_\varphi / 2` ordinates on the
        :math:`\xi > 0` half — each level a strictly η-monotone ARC in
        march order (:eq:`folded-level-arc`).

        **This is a QUOTIENT, not a half-range restriction**: each
        orbit's representative carries the WHOLE orbit weight
        (:math:`2w` — the action is free), so total mass stays
        :math:`4\pi` bit-exactly and no consumer carries a
        compensating factor.  ξ-odd functions are *not in the
        quotient's space* — the honest phase space for 1-D
        azimuthally-symmetric transport, where :math:`\psi` is ξ-even
        by geometric theorem.

        **The offset is DERIVED, not chosen** (T25): the fold is free
        — no node fixed by the mirror, :math:`\Sigma = \emptyset` —
        exactly when the azimuthal nodes straddle the mirror planes,
        i.e. the STAGGERED shift at even :math:`n_\varphi`.  That is
        why this factory exists rather than a ``fold=`` flag on
        :meth:`product`: offset-selection and fold-intent must be
        co-selected, or the incoherent NODE_ALIGNED-plus-fold
        combination (a node ON :math:`\Sigma`, an edge-node march
        start, the :math:`\tau = 0` division) becomes spellable.

        Consequences on the M-M closure (all gated): every level
        consumes an independent ψ½ seed
        (:func:`~orpheus.sn.sweep.pole_angular_closure.march_start_structure_per_level`
        — carrying), :math:`\tau_{\rm raw} \in [\tfrac15, \tfrac45]`
        strictly away from the :math:`\{0, 1\}` singularities, and the
        reversal identity :math:`\tau_m + \tau_{M-1-m} = 1` holds
        bit-exactly (:eq:`morel-montry-folded-arc`).

        Raises
        ------
        ValueError
            If ``n_phi`` is odd: the staggered rule at odd
            :math:`n_\varphi` places one node per level ON the mirror
            (:math:`\varphi = \pi \in \Sigma`), so the fold is not
            free and the arc would start on an edge node.
        """
        if n_phi % 2 != 0:
            raise ValueError(
                f"folded_product requires an even n_phi: the staggered "
                f"circle rule at odd n_phi={n_phi} places one node per "
                f"level ON the mirror plane (phi = pi is in Sigma), so "
                f"the sigma_y fold is not free and the folded arc would "
                f"start on an edge node (tau_raw = 0, the division the "
                f"fold exists to remove)."
            )
        measure, structure = spherical_product(
            gauss_legendre_on_mu(n_mu),
            periodic_trapezoid(n_phi, shift=STAGGERED),
        )
        full = cls(measure=measure, level_structure=structure)
        return full.quotient(SubgroupOfO3.Mirror("y"))


__all__ = ["Quadrature"]
