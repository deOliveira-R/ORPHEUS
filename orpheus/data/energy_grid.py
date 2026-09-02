r"""Energy-group structure — the energy-axis value object + its dual frame views.

:class:`EnergyGrid` is a multigroup energy structure: a strictly **DESCENDING**
boundary array (the canonical fast-first convention, group ``0`` = fastest; see
:ref:`canonical-group-convention`). ``N+1`` boundaries → ``N`` groups. It is the
energy analogue of a coarse :class:`~orpheus.geometry.mesh.Mesh1D`, and — exactly like
``Mesh1D`` — it yields **both** halves of a discrete frame:

* :meth:`EnergyGrid.as_measure` → a :class:`~orpheus.numerics.measure.DiscreteMeasure`,
  the **source** view (the axis you project FROM), symmetric with ``Mesh1D.volume_measure``;
* :meth:`EnergyGrid.as_basis` → an :class:`~orpheus.numerics.basis.IndicatorBasis`, the
  **target** view (the axis you project TO), symmetric with ``Mesh1D.indicator_basis``.

The group-structure **mismatch treatment** — when the source and target grids do not
align — is the *binary* :meth:`EnergyGrid.overlap_to` (it reads BOTH grids' edges, so it
cannot be a unary view): ``coarse.as_basis()`` (nested, one-hot) ⊂
``fine.overlap_to(coarse)`` (non-nested, fractional). The collapse VERB that consumes
these views is :meth:`orpheus.data.macro_xs.mixture.Mixture.condense`; this module owns
only the grid + its frame views (data-native — no transport dependency).

Why fractional overlap (the non-nested problem)
===============================================

The production case (a 421-group library → the WIMS-D 69/172 structures) is
**non-nested**: the coarse boundaries do not align with the fine grid, and the
coarse grid is even *locally finer* than the fine spacing in narrow resonance /
thermal bands. A naive "assign each fine group wholly to one coarse group" leaves
some coarse groups empty. The standard reactor-physics fix is fractional re-binning
(Hébert, *Applied Reactor Physics* §3.5; the flux-weighted average is the rank-0
moment of Generalized Energy Condensation, Rahnema et al. NSE 160:41, 2008): a
straddling fine group's contribution is split by the fraction of its (within-group
flux–weighted) interval falling in each coarse group,

.. math::

   f_{g,G} \;=\; \frac{\int_{g\cap G} w(E)\,\mathrm{d}E}{\int_g w(E)\,\mathrm{d}E},

with :math:`w(E)` a **selectable** within-group flux model
(:class:`WithinGroupSpectrum`; 1/E is the default, :class:`InverseEnergySpectrum`).
:meth:`overlap_to` returns the partition-of-unity table :math:`T[g,G]\in[0,1]` (rows
summing to 1) carried by an :class:`~orpheus.numerics.basis.OverlapBasis`; a *nested*
coarse grid gives the one-hot degenerate (the exact group-sum).

Downsampling only
-----------------

Condensation **loses** information; it is well-posed fine→coarse only. Reconstructing
a *finer* structure from group-integrated data fabricates sub-group detail the data
never carried, so :meth:`overlap_to` **refuses** a coarse target with more groups than
the source (the upscaling guard). Where the coarse grid is *locally* finer than the fine
grid, the within-group model performs a bounded, explicit interpolation — those coarse
groups are reported by :attr:`~orpheus.numerics.basis.OverlapBasis.fractional_columns`
so the data-vs-assumption provenance is visible. (Faithful reconstruction / honest
upscaling via rank>0 GEC moments is a future capability — GitHub #275.)

The reuse — index-space membership, the existing PG frame
=========================================================

The fractional table feeds the *unchanged* Petrov-Galerkin frame
(:class:`~orpheus.numerics.frame.PetrovGalerkinFrame`,
:class:`~orpheus.numerics.basis.OverlapBasis`,
:class:`~orpheus.numerics.basis.WeightedIndicatorBasis`): the flux is the **test
weight**, the energy measure is **counting** (``w_g = 1`` — :meth:`as_measure`; ORPHEUS
``φ_g`` is already group-integrated). Because the table is a partition of unity
(:attr:`~orpheus.numerics.basis.base.GramStructure.PARTITION_OF_UNITY`),
``frame.project = G⁻¹ M`` yields the rate-preserving average with the diagonal Gram
:math:`\Phi_G = \sum_g T[g,G]\,\varphi_g` — no new frame method, only the
fractional-overlap basis.

References
----------

* Hébert, A. (2009). *Applied Reactor Physics*. Polytechnique. §3.5 (Eq. 3.103–3.112
  — multigroup condensation: vector flux-weighted, scattering 2-axis, χ pure sum).
* Rahnema, F., Douglass, S., Forget, B. (2008). Generalized Energy Condensation
  Theory. *Nucl. Sci. Eng.* 160:41. DOI 10.13182/NSE160-41 (the flux-weighted
  average is the rank-0 moment).
* IAEA-TECDOC / STI/Pub/1264 (2007). *WIMS-D Library Update*, Tables 11.1–11.3.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from typing import Protocol, runtime_checkable

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.basis.indicator_basis import IndicatorBasis
from orpheus.numerics.basis.overlap_basis import OverlapBasis
from orpheus.numerics.manifold import EnergyGroups
from orpheus.numerics.measure import DiscreteMeasure

__all__ = [
    "EnergyGrid",
    "InverseEnergySpectrum",
    "WithinGroupSpectrum",
]

# Relative tolerance for the coarse-within-fine span check (FP slack on a boundary
# that should coincide with the fine ceiling).
_SPAN_RTOL = 1e-9


@dataclass(frozen=True, eq=False)
class EnergyGrid:
    r"""A multigroup energy structure: strictly DESCENDING boundary energies.

    The canonical ORPHEUS energy axis (see :ref:`canonical-group-convention`): group
    index ``0`` is the FASTEST (highest energy), so the boundary array is strictly
    decreasing — ``edges[0]`` the highest bound, ``edges[-1]`` the lowest (the floor,
    which may be ``0`` for a thermal structure such as WIMS-69). ``N+1`` boundaries
    define ``N`` groups; group ``g`` occupies ``[edges[g+1], edges[g])``.

    Parameters
    ----------
    edges : NDArray
        The ``(n_groups + 1,)`` boundary energies in eV, strictly DESCENDING,
        non-negative. Coerced to a float array at construction.

    Notes
    -----
    Frozen, ``eq=False`` (identity equality; the sole field is a NumPy array). The
    strict-monotonicity + positivity invariant is the #265 data-layer value-object
    slice — checked once at construction, never re-checked (frozen).
    """

    edges: NDArray

    def __post_init__(self) -> None:
        edges = np.ascontiguousarray(self.edges, dtype=float)
        object.__setattr__(self, "edges", edges)
        if edges.ndim != 1 or edges.shape[0] < 2:
            raise ValueError(
                f"EnergyGrid needs a 1-D array of ≥2 boundaries, got shape "
                f"{edges.shape}"
            )
        if not bool(np.all(edges[:-1] > edges[1:])):
            raise ValueError(
                "EnergyGrid boundaries must be strictly DESCENDING (group 0 = "
                "fastest, the canonical fast-first convention); got "
                f"edges={edges!r}"
            )
        if edges[-1] < 0.0:
            raise ValueError(
                f"energy boundaries must be non-negative (a thermal floor of 0 eV "
                f"is allowed); got floor edges[-1]={edges[-1]!r}"
            )

    @property
    def n_groups(self) -> int:
        """The group count :math:`N` (``len(edges) - 1``)."""
        return int(self.edges.shape[0] - 1)

    @cached_property
    def representative_energy(self) -> NDArray:
        r"""Per-group representative energy — the geometric midpoint (descending).

        :math:`\bar E_g = \sqrt{E_g^{\mathrm{upper}}\, E_g^{\mathrm{lower}}}`, the
        natural centre on the log-spaced energy axis. The bottom group's lower edge
        may be ``0`` (a thermal floor); there the geometric mean degenerates, so the
        representative falls back to half the upper edge (still strictly inside the
        group).
        """
        upper = self.edges[:-1]
        lower = self.edges[1:]
        return np.where(lower > 0.0, np.sqrt(upper * lower), 0.5 * upper)

    @cached_property
    def energy_widths(self) -> NDArray:
        r"""Per-group energy width :math:`\Delta E_g = E_g^{\mathrm{upper}} - E_g^{\mathrm{lower}}`.

        Positive by construction (edges are descending). This is the denominator
        of the **flux-per-unit-energy** spectrum :math:`\phi_g / \Delta E_g` — the
        single source of the group widths the homogeneous / MoC diagnostics plot.
        """
        return self.edges[:-1] - self.edges[1:]

    @cached_property
    def lethargy_widths(self) -> NDArray:
        r"""Per-group lethargy width :math:`\Delta u_g = \ln\!\bigl(E_g^{\mathrm{upper}} / E_g^{\mathrm{lower}}\bigr)`.

        Positive by construction (edges are descending). The denominator of the
        **flux-per-unit-lethargy** spectrum :math:`\phi_g / \Delta u_g`. A thermal
        floor (``edges[-1] == 0``) gives :math:`+\infty` for the bottom group — its
        lethargy width is genuinely unbounded (:math:`u \to \infty` as
        :math:`E \to 0`), unlike the geometric :attr:`representative_energy`, which
        degenerates there and falls back to half the upper edge.
        """
        return np.log(self.edges[:-1] / self.edges[1:])

    # ── The dual frame views (source = measure, target = basis) ───────────
    def as_measure(self) -> DiscreteMeasure:
        r"""The energy-axis **counting** measure — the SOURCE view (project FROM).

        Nodes are the integer group indices ``0…n_groups-1``; weights are all ``1``
        because ORPHEUS ``φ_g`` is already group-integrated, so the discrete reaction
        rate is a plain group sum (a lethargy/Δu weight would break rate preservation —
        the within-group ``w`` lives in the overlap fractions, not here).
        ``support=EnergyGroups()`` IS the physical phase-space factor
        (:attr:`~orpheus.numerics.measure.DiscreteMeasure.phase` reads it). Symmetric
        with :meth:`~orpheus.geometry.mesh.Mesh1D.volume_measure` on the spatial axis.
        """
        n = self.n_groups
        return DiscreteMeasure(
            # ⭐ The SAME manifold :meth:`as_basis` partitions, three lines
            # below — one axis, one point set. The two used to say it
            # separately (``EnergyGroups()`` here, the string ``"energy"``
            # there), which is the divergence ``test_d6`` exists to refuse.
            nodes=np.arange(n, dtype=float),
            weights=np.ones(n),
            # ⚠ ``EnergyGroups(n)``, not the bare ``EnergyGroups()``: the group
            # COUNT is part of the point set's identity, and this measure knows
            # it. Under the retired string tag both halves said ``"energy"`` and
            # agreed by being equally uninformative — the tag was lossy, and the
            # loss is what hid the disagreement (``as_basis`` has always built
            # ``EnergyGroups(ng=n)``).
            support=EnergyGroups(n),
        )

    def as_basis(self) -> IndicatorBasis:
        r"""The group-indicator basis in index space — the TARGET view (project TO).

        A plain :class:`~orpheus.numerics.basis.IndicatorBasis` over the group-index
        partition (edges ``-0.5, 0.5, …, n-0.5``) — the **nested/one-hot** target view.
        The non-nested fractional target is :meth:`overlap_to` (an
        :class:`~orpheus.numerics.basis.OverlapBasis`, which IS-A ``IndicatorBasis``
        carrying the mismatch table). Symmetric with
        :meth:`~orpheus.geometry.mesh.Mesh1D.indicator_basis` on the spatial axis.
        """
        index_edges = np.arange(self.n_groups + 1, dtype=float) - 0.5
        return IndicatorBasis(
            edges_per_axis=(index_edges,),
            partition_of=EnergyGroups(self.n_groups),
        )

    def overlap_to(
        self, coarse: "EnergyGrid", /, within_group: "WithinGroupSpectrum | None" = None,
    ) -> OverlapBasis:
        r"""The fine→coarse fractional-overlap trial basis — the BINARY mismatch treatment.

        ``self`` is the fine (source) grid, ``coarse`` the target. Returns the
        :class:`~orpheus.numerics.basis.OverlapBasis` carrying the partition-of-unity
        table :math:`T[g,G]\in[0,1]`: a fine group straddling a coarse boundary
        apportions a fraction of its rate to each coarse group it overlaps, weighted by
        ``within_group`` (1/E by default, :class:`InverseEnergySpectrum`). A *nested*
        coarse grid gives the one-hot degenerate (``coarse.as_basis()`` carrying an
        identity table).

        This is the irreducibly **binary** ``(fine, coarse)`` step — it reads BOTH grids'
        edges, so it cannot be a unary view: ``coarse.as_basis()`` (nested, one-hot) ⊂
        ``fine.overlap_to(coarse)`` (non-nested, fractional). It is the data-level home
        the retired ``GroupCondensation`` collapsed into.

        Raises
        ------
        ValueError
            If ``coarse``'s ceiling exceeds ``self``'s (the coarse structure claims
            energy coverage the fine data lacks), or if ``coarse`` has MORE groups than
            ``self`` (the upscaling guard — condensation only downsamples; faithful
            reconstruction is GitHub #275).
        """
        w = within_group if within_group is not None else InverseEnergySpectrum()
        if coarse.edges[0] > self.edges[0] * (1.0 + _SPAN_RTOL):
            raise ValueError(
                f"coarse ceiling {coarse.edges[0]:.4e} eV exceeds the fine ceiling "
                f"{self.edges[0]:.4e} eV — ill-posed: the coarse structure claims energy "
                f"coverage the fine data does not have."
            )
        if coarse.n_groups > self.n_groups:
            raise ValueError(
                f"upscaling refused: coarse has {coarse.n_groups} groups > fine "
                f"{self.n_groups}. Condensation only DOWNSAMPLES; a finer target would "
                f"fabricate sub-group structure the data does not contain (faithful "
                f"reconstruction is GitHub #275)."
            )
        table = self._overlap_table(coarse, w)
        if not bool(np.allclose(table.sum(axis=1), 1.0, atol=1e-9)):
            raise ValueError(
                "internal error: the fractional-overlap table is not a partition of "
                "unity (every fine group's weight must sum to 1 across coarse groups)."
            )
        # The trial IS the coarse grid's basis-view (the target you project ONTO),
        # decorated with the fractional overlap table (the mismatch treatment).
        # The trial EATS this (fine) grid's groups and SPANS the coarse ones.
        return OverlapBasis.from_indicator(
            coarse.as_basis(), table, fine=EnergyGroups(self.n_groups),
        )

    def _overlap_table(
        self, coarse: "EnergyGrid", within_group: "WithinGroupSpectrum", /,
    ) -> NDArray:
        r"""The ``(n_fine, n_coarse)`` fractional membership table :math:`T[g,G]\in[0,1]`.

        ``T[g,G]`` is the fraction of fine group ``g``'s (``within_group``-weighted)
        interval lying in coarse group ``G``; rows sum to 1. The outer coarse edges are
        clamped to the fine span so every fine group's full weight is captured
        (over-ceiling fine mass → coarse 0; under-floor → last).
        """
        fe = self.edges
        n_fine = self.n_groups
        n_coarse = coarse.n_groups
        clamped_edges = coarse.edges.astype(float).copy()
        clamped_edges[0] = max(clamped_edges[0], fe[0])
        clamped_edges[-1] = min(clamped_edges[-1], fe[-1])
        coarse_hi = clamped_edges[:-1]
        coarse_lo = clamped_edges[1:]
        table = np.zeros((n_fine, n_coarse))
        for g in range(n_fine):
            g_hi = fe[g]
            g_lo = fe[g + 1]
            lo = np.maximum(g_lo, coarse_lo)
            hi = np.minimum(g_hi, coarse_hi)
            mask = hi > lo
            if bool(np.any(mask)):
                table[g, mask] = within_group.integrated_weight(
                    lo[mask], hi[mask]
                ) / within_group.integrated_weight(np.asarray(g_lo), np.asarray(g_hi))
        return table


# ──────────────────────────────────────────────────────────────────────────────
# Within-group flux model (the selectable straddle-apportionment strategy)
# ──────────────────────────────────────────────────────────────────────────────


@runtime_checkable
class WithinGroupSpectrum(Protocol):
    r"""The assumed sub-fine-group flux shape :math:`w(E)`.

    Apportioning a fine group that straddles a coarse boundary needs a within-group
    flux model: the fraction going to a coarse group is
    :math:`f_{g,G} = \big(\int_{g\cap G} w\big) / \big(\int_g w\big)`. Implementations
    provide the definite integral over an energy interval; the fraction (and hence
    the partition of unity) follows by the integral's additivity.
    """

    def integrated_weight(self, lo: NDArray, hi: NDArray) -> NDArray:
        r""":math:`\int_{lo}^{hi} w(E)\,\mathrm{d}E`, vectorised over interval arrays."""
        ...


@dataclass(frozen=True)
class InverseEnergySpectrum:
    r"""Within-group flux model :math:`w(E)\propto 1/E` — flat in lethargy (the default).

    The asymptotic slowing-down spectrum: :math:`\int_{lo}^{hi} \mathrm{d}E/E =
    \ln(hi/lo)` is the **lethargy width**, so a straddling fine group's overlap
    fraction is a lethargy ratio. The standard first choice for condensation
    (NJOY ``IWT=3``). Requires positive energies (1/E diverges at ``E=0``); a
    zero-floor *source* grid needs a thermal-cutoff or a different model.
    """

    def integrated_weight(self, lo: NDArray, hi: NDArray) -> NDArray:
        lo_arr = np.asarray(lo, dtype=float)
        hi_arr = np.asarray(hi, dtype=float)
        if bool(np.any(lo_arr <= 0.0)):
            raise ValueError(
                "InverseEnergySpectrum (1/E) is undefined at E≤0; the fine grid "
                "must have positive group boundaries (a zero-floor grid needs a "
                "thermal-cutoff or a different WithinGroupSpectrum)."
            )
        return np.log(hi_arr / lo_arr)
