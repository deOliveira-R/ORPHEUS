r"""Energy-group structure and condensation partition — the energy-axis value objects.

Energy condensation collapses fine-group cross sections onto a coarse group
structure, spectrum-weighted (the energy-axis transpose of spatial homogenisation;
see :meth:`orpheus.sn.solution.Solution.condense`). Two value objects carry its
structure:

* :class:`EnergyGrid` — a multigroup energy structure: a strictly **DESCENDING**
  boundary array (the canonical fast-first convention, group ``0`` = fastest; see
  :ref:`canonical-group-convention`). ``N+1`` boundaries → ``N`` groups. The energy
  analogue of a coarse :class:`~orpheus.geometry.mesh.Mesh1D`.
* :class:`GroupCondensation` — the fine→coarse partition map. Built by
  **conservative fractional re-binning**: a fine group straddling a coarse boundary
  apportions a *fraction* of its rate to each coarse group it overlaps, so the
  membership table :math:`T[g,G]\in[0,1]` is a **partition of unity** (rows sum to
  1). A *nested* coarse grid (boundaries a subset of the fine grid) gives the one-hot
  degenerate and the exact group-sum.

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

Downsampling only
-----------------

Condensation **loses** information; it is well-posed fine→coarse only. Reconstructing
a *finer* structure from group-integrated data fabricates sub-group detail the data
never carried, so :class:`GroupCondensation` **refuses** a coarse target with more
groups than the source (the upscaling guard). Where the coarse grid is *locally*
finer than the fine grid, the within-group model performs a bounded, explicit
interpolation — those coarse groups are reported by
:attr:`~GroupCondensation.locally_interpolated` so the data-vs-assumption provenance
is visible. (Faithful reconstruction / honest upscaling via rank>0 GEC moments is a
future capability — GitHub #275.)

The reuse — index-space membership, the existing PG frame
=========================================================

The fractional table feeds the *unchanged* Petrov-Galerkin frame
(:class:`~orpheus.numerics.frame.PetrovGalerkinFrame`,
:class:`~orpheus.numerics.basis.OverlapBasis`,
:class:`~orpheus.numerics.basis.WeightedIndicatorBasis`): the flux is the **test
weight**, the energy measure is **counting** (``w_g = 1``; ORPHEUS ``φ_g`` is already
group-integrated). Because the table is a partition of unity,
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

from dataclasses import dataclass, field
from functools import cached_property
from typing import Protocol, runtime_checkable

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.basis.overlap_basis import OverlapBasis
from orpheus.numerics.measure import DiscreteMeasure

__all__ = [
    "EnergyGrid",
    "GroupCondensation",
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

    def condense_to(
        self, coarse: "EnergyGrid", /, within_group: "WithinGroupSpectrum | None" = None
    ) -> "GroupCondensation":
        """Derive the fine→coarse :class:`GroupCondensation` onto ``coarse``.

        ``self`` is the fine grid. ``within_group`` selects the sub-fine-group flux
        model for straddle apportionment (default 1/E, :class:`InverseEnergySpectrum`).
        """
        if within_group is None:
            return GroupCondensation(self, coarse)
        return GroupCondensation(self, coarse, within_group)


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


# ──────────────────────────────────────────────────────────────────────────────
# The fine→coarse partition (fractional overlap)
# ──────────────────────────────────────────────────────────────────────────────


@dataclass(frozen=True, eq=False)
class GroupCondensation:
    r"""The fine→coarse group partition: a fractional, partition-of-unity overlap map.

    Built by conservative fractional re-binning — a fine group straddling a coarse
    boundary apportions a fraction of its rate to each coarse group it overlaps
    (weighted by ``within_group``). The membership :attr:`table` is a partition of
    unity (rows sum to 1); a *nested* coarse grid gives the one-hot degenerate and
    the exact group-sum.

    Parameters
    ----------
    fine, coarse : EnergyGrid
        The fine (source) and coarse (target) energy structures.
    within_group : WithinGroupSpectrum, optional
        The sub-fine-group flux model for straddle apportionment. Default 1/E
        (:class:`InverseEnergySpectrum`); selectable.

    Raises
    ------
    ValueError
        If ``coarse``'s ceiling exceeds ``fine``'s (the coarse structure claims
        coverage the data lacks), or if ``coarse`` has MORE groups than ``fine`` (the
        upscaling guard — condensation only downsamples).

    Notes
    -----
    Frozen, ``eq=False``. Over-ceiling / under-floor fine mass clamps into the first /
    last coarse group, so the partition of unity holds for every fine group.
    """

    fine: EnergyGrid
    coarse: EnergyGrid
    within_group: WithinGroupSpectrum = field(default_factory=InverseEnergySpectrum)

    def __post_init__(self) -> None:
        if self.coarse.edges[0] > self.fine.edges[0] * (1.0 + _SPAN_RTOL):
            raise ValueError(
                f"coarse ceiling {self.coarse.edges[0]:.4e} eV exceeds the fine "
                f"ceiling {self.fine.edges[0]:.4e} eV — ill-posed: the coarse "
                f"structure claims energy coverage the fine data does not have."
            )
        if self.coarse.n_groups > self.fine.n_groups:
            raise ValueError(
                f"upscaling refused: coarse has {self.coarse.n_groups} groups > fine "
                f"{self.fine.n_groups}. Condensation only DOWNSAMPLES; a finer target "
                f"would fabricate sub-group structure the data does not contain "
                f"(faithful reconstruction is GitHub #275)."
            )
        row_sums = self.table.sum(axis=1)  # triggers table + within-group validation
        if not bool(np.allclose(row_sums, 1.0, atol=1e-9)):
            raise ValueError(
                "internal error: the fractional-overlap table is not a partition of "
                "unity (every fine group's weight must sum to 1 across coarse groups)."
            )

    @cached_property
    def table(self) -> NDArray:
        r"""The fine→coarse fractional membership table :math:`T[g,G]\in[0,1]`.

        ``(n_fine, n_coarse)``, rows summing to 1. ``T[g,G]`` is the fraction of fine
        group ``g``'s (``within_group``-weighted) interval lying in coarse group
        ``G``. The :math:`@T` contraction is the **sink-sum** half of the scattering
        2-axis collapse and the pure birth-group sum for ``χ`` (``χ_G = χ @ T``);
        ``frame.project`` against it does the source-axis flux-average. Equal to
        ``self.indicator_basis().evaluate(self.measure.nodes)`` by construction.
        """
        fe = self.fine.edges
        n_fine = self.fine.n_groups
        n_coarse = self.coarse.n_groups
        # Clamp the OUTER coarse edges to the fine span so every fine group's full
        # weight is captured (over-ceiling fine mass → coarse 0; under-floor → last).
        ce = self.coarse.edges.astype(float).copy()
        ce[0] = max(ce[0], fe[0])
        ce[-1] = min(ce[-1], fe[-1])
        coarse_hi = ce[:-1]
        coarse_lo = ce[1:]
        w = self.within_group
        table = np.zeros((n_fine, n_coarse))
        for g in range(n_fine):
            g_hi = fe[g]
            g_lo = fe[g + 1]
            lo = np.maximum(g_lo, coarse_lo)
            hi = np.minimum(g_hi, coarse_hi)
            mask = hi > lo
            if bool(np.any(mask)):
                table[g, mask] = w.integrated_weight(
                    lo[mask], hi[mask]
                ) / w.integrated_weight(np.asarray(g_lo), np.asarray(g_hi))
        return table

    @cached_property
    def coarse_of_fine(self) -> NDArray:
        r"""The DOMINANT coarse group of each fine group (``argmax`` of :attr:`table`).

        For a nested condensation the table is one-hot, so this is the exact
        containing-coarse-group map (and the partition is contiguous, fast-first).
        For a straddling fine group it is the coarse group receiving the largest
        fraction — a reporting/diagnostic view; the real apportionment is
        :attr:`table`.
        """
        return self.table.argmax(axis=1).astype(int)

    @cached_property
    def locally_interpolated(self) -> NDArray:
        r"""Coarse-group indices whose data leans on the within-group flux model.

        The coarse groups that received a **fractional** (strictly between 0 and 1)
        contribution from any fine group — i.e. where a fine group straddled the
        boundary and ``within_group`` performed local interpolation. Empty for a
        nested condensation (pure rate-preserving collapse, no assumption); non-empty
        where the coarse grid is locally finer than the fine grid. The provenance
        report of the downsampling/interpolation asymmetry.
        """
        frac = (self.table > 1e-12) & (self.table < 1.0 - 1e-12)
        return np.nonzero(frac.any(axis=0))[0]

    def indicator_basis(self) -> OverlapBasis:
        r"""The coarse trial basis carrying the fractional :attr:`table`.

        An :class:`~orpheus.numerics.basis.OverlapBasis` (the partition-of-unity
        generalisation of :class:`~orpheus.numerics.basis.IndicatorBasis`); its
        ``edges_per_axis`` is the coarse-group index partition (for the coarse cell
        count + coefficient space) and its membership is the precomputed overlap
        table. Feeds the unchanged ``frame.project``.
        """
        n_coarse = self.coarse.n_groups
        coarse_index_edges = np.arange(n_coarse + 1, dtype=float) - 0.5
        return OverlapBasis(
            edges_per_axis=(coarse_index_edges,), overlap_table=self.table
        )

    @property
    def measure(self) -> DiscreteMeasure:
        r"""The fine energy-axis **counting** measure (``weights = 1``).

        Nodes are the integer fine-group indices ``0…ng-1`` (the bucketing coordinate
        for :meth:`indicator_basis`); weights are all ``1`` because ORPHEUS ``φ_g`` is
        already group-integrated, so the discrete reaction rate is a plain group sum
        (a lethargy/Δu weight would break rate preservation — the within-group ``w``
        lives in the overlap fractions, not here). ``support="energy"`` tags the
        physical phase-space factor (:attr:`DiscreteMeasure.phase` reads it).
        """
        n_fine = self.fine.n_groups
        return DiscreteMeasure(
            nodes=np.arange(n_fine, dtype=float),
            weights=np.ones(n_fine),
            support="energy",
        )

    @classmethod
    def from_grids(
        cls,
        fine: EnergyGrid,
        coarse: EnergyGrid,
        /,
        within_group: "WithinGroupSpectrum | None" = None,
    ) -> "GroupCondensation":
        """Build the condensation from fine and coarse :class:`EnergyGrid` instances.

        The canonical call-site spelling; the fractional-overlap derivation +
        well-posedness checks run in ``__post_init__``.
        """
        if within_group is None:
            return cls(fine, coarse)
        return cls(fine, coarse, within_group)
