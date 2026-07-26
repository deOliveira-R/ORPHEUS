"""Macroscopic cross section computation for isotope mixtures.

A Mixture holds the isotopes, their number densities, and the resulting
macroscopic cross sections used by reactor solvers.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Optional

import numpy as np
from scipy.sparse import csr_matrix

from orpheus.data.emission_spectrum import enforce_emission_spectrum
from orpheus.data.micro_xs.isotope import NG, Isotope
from .interpolation import interp_sig_s, interp_xs_field
from .sigma_zeros import solve_sigma_zeros

if TYPE_CHECKING:
    from orpheus.data.energy_grid import EnergyGrid, WithinGroupSpectrum


@dataclass
class Mixture:
    """Macroscopic cross sections for a homogeneous mixture.

    Attributes
    ----------
    SigC, SigL, SigF, SigP, SigT : (NG,) arrays
        Macroscopic capture, (n,alpha), fission, production, total XS in 1/cm.
    SigS : list of (NG, NG) sparse matrices, one per Legendre order.
    Sig2 : (NG, NG) sparse — macroscopic (n,2n) matrix.
    chi  : (NG,) — fission spectrum of the mixture.
    eg   : (NG+1,) energy group boundaries in eV, *or* ``None``.

        For real production cases (XS computed via :func:`compute_macro_xs`
        from ENDF :class:`~orpheus.data.micro_xs.isotope.Isotope`), ``eg``
        is populated from ``isotopes[0].eg`` — the energy grid is part of
        the upstream nuclear-data library and is genuinely defined.

        For synthetic / Sood-style abstract cross sections (built via
        :func:`orpheus.derivations.common.xs_library.make_mixture` or
        constructed directly in MMS / verification harnesses), ``eg`` is
        ``None``: there is no real grid. Downstream consumers (lethargy
        widths, per-energy plotting, spectrum-weighted condensation) MUST
        guard for ``None`` and either skip the per-energy path or raise
        a clear error. ``None`` is the intended representation when no
        physical energy grid exists — *not* a placeholder to be filled in.

        ``eg``'s production use is spectrum-weighted energy condensation
        (Petrov-Galerkin fractional-overlap projection from fine to coarse
        groups; see :meth:`condense`). That step depends on the grid; it
        cannot run on synthetic mixtures.
    """

    SigC: np.ndarray
    SigL: np.ndarray
    SigF: np.ndarray
    SigP: np.ndarray
    SigT: np.ndarray
    SigS: list[csr_matrix]
    Sig2: csr_matrix
    chi: np.ndarray
    eg: np.ndarray | None = None

    def __post_init__(self) -> None:
        # Coerce chi to the validated value-object and enforce the simplex
        # / null law at the data source (mirrors Isotope.__post_init__). χ
        # is consumed only as a fission SOURCE (χ·νΣ_f·φ), so the law keys
        # on PRODUCTION (SigP = νΣ_f > 0): a producing mixture's spectrum is
        # a probability simplex; a non-producing mixture emits no fission
        # neutrons, so its spectrum is null.
        self.chi = enforce_emission_spectrum(self.chi, is_producing=self.is_producing)

    @property
    def is_producing(self) -> bool:
        r"""``True`` iff the mixture emits fission neutrons (:math:`\nu\Sigma_f > 0`).

        Keys on the production XS :attr:`SigP` (= :math:`\nu\Sigma_f`), the
        quantity that actually drives the fission source ``χ·SigP·φ``. This
        is the predicate the χ emission-spectrum law keys on.
        """
        return bool(np.any(self.SigP > 0))

    @property
    def ng(self) -> int:
        """Number of energy groups (inferred from data)."""
        return len(self.SigT)

    @property
    def energy_grid(self) -> "EnergyGrid":
        r"""This mixture's source :class:`~orpheus.data.energy_grid.EnergyGrid` (from :attr:`eg`).

        The XS carries its own energy structure, so it carries a measure (the SOURCE
        view a condensation projects FROM — :meth:`EnergyGrid.as_measure`). Raises if
        ``eg is None`` (a synthetic/Sood mixture has no physical grid to condense).
        """
        if self.eg is None:
            raise ValueError(
                "Mixture has no energy grid (eg is None — a synthetic/Sood mixture has "
                "no physical grid); cannot view it as an EnergyGrid."
            )
        from orpheus.data.energy_grid import EnergyGrid

        return EnergyGrid(self.eg)

    @property
    def absorption_xs(self) -> np.ndarray:
        """(NG,) absorption XS: fission + capture + (n,alpha) + (n,2n) out."""
        return self.SigF + self.SigC + self.SigL + np.array(self.Sig2.sum(axis=1)).ravel()

    @property
    def out_scattering_xs(self) -> np.ndarray:
        """(NG,) out-of-group scattering XS (P0 off-diagonal sum)."""
        return self.total_scattering_xs - self.in_scattering_xs

    @property
    def in_scattering_xs(self) -> np.ndarray:
        """(NG,) in-group (elastic) scattering XS (P0 diagonal)."""
        return np.array(self.SigS[0].diagonal()).ravel()

    @property
    def total_scattering_xs(self) -> np.ndarray:
        """(NG,) total scattering XS (P0 row sum = in + out)."""
        return np.array(self.SigS[0].sum(axis=1)).ravel()

    @property
    def n2n_out_xs(self) -> np.ndarray:
        """(NG,) total (n,2n) out-scattering XS (Sig2 row sum).

        Mirrors :attr:`total_scattering_xs` for the (n,2n) channel so the
        balance identity reads like the conservation law it encodes (see
        :meth:`balance_residual`).
        """
        return np.array(self.Sig2.sum(axis=1)).ravel()

    @property
    def transport_xs(self) -> np.ndarray:
        r"""(NG,) transport cross section (outflow approximation) in 1/cm.

        .. math::

           \Sigma_{tr,g} \;=\; \Sigma_{t,g} \;-\; \sum_{g'} \Sigma_{s1,\,g\to g'}

        — the total XS minus the P1 out-scatter row sum: the *outflow*
        (out-scatter) transport approximation (Stamm'ler & Abbate 1983).
        For a single medium with :math:`\Sigma_{s1} = \bar\mu\,\Sigma_{s0}`
        this is the textbook :math:`\Sigma_{tr} = \Sigma_t - \bar\mu\,
        \Sigma_s` (Duderstadt & Hamilton 1976). The flux-weighted *inflow*
        (in-scatter) refinement needs a solved spectrum and is therefore
        deliberately NOT a static data property.

        When the mixture carries no P1 moment (``len(SigS) == 1`` — the
        synthetic / Sood-style P0-only path), the out-scatter row sum is
        identically zero and :math:`\Sigma_{tr} = \Sigma_t` EXACTLY: the
        correct isotropic-scattering limit, not a fallback.

        Consumers: :attr:`diffusion_coefficient` (the diffusion data seam,
        #290); any future transport-corrected P0 treatment (CP/MoC).
        """
        if len(self.SigS) > 1:
            p1_out = np.array(self.SigS[1].sum(axis=1)).ravel()
        else:
            p1_out = np.zeros(self.ng)
        return self.SigT - p1_out

    @property
    def diffusion_coefficient(self) -> np.ndarray:
        r"""(NG,) diffusion coefficient :math:`D_g = 1/(3\,\Sigma_{tr,g})` in cm.

        The Fick's-law coefficient of the multigroup diffusion model
        (:eq:`diffusion-coefficient`), built on the outflow
        :attr:`transport_xs`.

        Raises
        ------
        ValueError
            If any group's :math:`\Sigma_{tr} \le 0` — a P1 out-scatter
            exceeding the total XS is unphysical data, and the reciprocal
            would silently produce a negative diffusion coefficient.
        """
        sig_tr = self.transport_xs
        if np.any(sig_tr <= 0.0):
            bad_group = int(np.argmin(sig_tr))
            raise ValueError(
                f"diffusion_coefficient requires transport_xs > 0 in every "
                f"group; got transport_xs[{bad_group}] = {sig_tr[bad_group]:g} "
                f"(P1 out-scatter row sum >= total XS is unphysical data)."
            )
        return 1.0 / (3.0 * sig_tr)

    @property
    def balance_residual(self) -> np.ndarray:
        r"""(NG,) per-group total-XS balance residual.

        The total cross section is, by definition, the sum of every removal
        channel:

        .. math::

           \Sigma_{t,g} = \Sigma_{c,g} + \Sigma_{L,g} + \Sigma_{f,g}
                          + \sum_{g'}\Sigma_{s0,g\to g'}
                          + \sum_{g'}\Sigma_{2n,g\to g'}

        i.e. capture + (n,alpha) + fission + total P0 scatter-out + (n,2n)
        out. This is VERBATIM the line that derives ``SigT`` in
        :func:`compute_macro_xs`. The residual is the per-group absolute
        defect :math:`|\Sigma_{t,g} - \text{(removal channels)}|`.

        ``SigP`` (production = :math:`\nu\Sigma_f`) is a fission-SOURCE
        multiplier, NOT a removal channel — it is deliberately absent.
        """
        derived = (
            self.SigC
            + self.SigL
            + self.SigF
            + self.total_scattering_xs
            + self.n2n_out_xs
        )
        return np.abs(self.SigT - derived)

    def assert_balanced(self, atol: float = 1e-9) -> None:
        r"""Assert the definitional total-XS identity holds per group.

        Raises :class:`ValueError` if the maximum per-group
        :attr:`balance_residual` exceeds ``atol``, reporting the residual
        and the offending group index.

        This is the canonical whole-identity check
        :math:`\Sigma_t = \Sigma_c + \Sigma_L + \Sigma_f
        + \text{rowsum}(\Sigma_{s0}) + \text{rowsum}(\Sigma_{2n})`.

        It is NOT enforced in ``__post_init__``: manufactured / synthetic
        cross sections legitimately carry non-physical totals — the Atalay
        1997 criticality encoding (``νΣ_f = (c-1)Σ_t`` with ``Σ_f = 0`` for
        ``c > 1``), structural test scaffolds, and the billiard ``SigP``
        carrier all build imbalanced ``Mixture`` instances on purpose. The
        law is invoked by PHYSICAL builders only (e.g. :func:`compute_macro_xs`,
        which derives ``SigT`` via the identity and therefore always balances —
        the call there is a free regression guard against a future derivation
        bug).
        """
        residual = self.balance_residual
        max_residual = float(residual.max())
        if max_residual > atol:
            bad_group = int(np.argmax(residual))
            raise ValueError(
                f"Mixture XS imbalance: max |SigT - (SigC+SigL+SigF"
                f"+rowsum(SigS0)+rowsum(Sig2))| = {max_residual:g} > atol={atol:g} "
                f"in group {bad_group}"
            )

    @property
    def scattering_ratio(self) -> np.ndarray:
        r"""(NG,) Case–Zweifel secondaries-per-collision parameter.

        Per-group :math:`c_g = (\Sigma_{s,g} + \nu\Sigma_{f,g}) / \Sigma_{t,g}`,
        where :math:`\Sigma_{s,g}` is the *total* P0 scattering out of
        group :math:`g` (in + out) and :math:`\nu\Sigma_{f,g}` is the
        production XS already stored as :attr:`SigP`. This is the
        Case–Zweifel ``c`` — the mean number of secondaries (scatter +
        fission neutrons) emitted per collision in group :math:`g`.
        For 1-group isotropic cases used by Sood/LA-13511 benchmarks,
        ``scattering_ratio[0]`` is the canonical :math:`c` parameter
        consumed by F_N and other analytical solvers.

        For purely scattering media this reduces to
        :math:`\Sigma_s/\Sigma_t \le 1`. For multiplying media
        (:math:`\nu\Sigma_f > 0`) the ratio can exceed 1 (the F_N
        c-sweep spans :math:`c \in [0.9, 1.1]` for sub-/super-critical
        configurations).

        Notes
        -----
        For multi-group, ``scattering_ratio[g]`` is a per-group
        diagnostic — solvers use the full (``SigS``, ``SigP``)
        operators rather than this scalar. The property exists for
        the 1G analytical bridge and for diagnostic display.
        """
        return (self.total_scattering_xs + self.SigP) / self.SigT

    @classmethod
    def from_dense_channels(
        cls,
        *,
        SigC: np.ndarray,
        SigL: np.ndarray,
        SigF: np.ndarray,
        SigP: np.ndarray,
        SigT: np.ndarray,
        SigS: list[np.ndarray],
        Sig2: np.ndarray,
        chi: np.ndarray,
        eg: np.ndarray | None,
    ) -> "Mixture":
        r"""Assemble a :class:`Mixture` from DENSE collapsed channels — the shared assembler.

        The single home for "build a Mixture from coarsened dense channel arrays": wraps
        the matrix channels (``SigS[ℓ]``, ``Sig2``) in :class:`~scipy.sparse.csr_matrix`
        and threads ``eg``. BOTH coarsening verbs — :meth:`condense` (energy) and
        :meth:`orpheus.transport.mesh.material_xs_field.MaterialXSField.project_through`
        (space) — call this, so the channel-assembly taxonomy lives once (Cardinal Rule 2)
        rather than duplicated per verb.
        """
        return cls(
            SigC=SigC, SigL=SigL, SigF=SigF, SigP=SigP, SigT=SigT,
            SigS=[csr_matrix(s) for s in SigS],
            Sig2=csr_matrix(Sig2),
            chi=chi, eg=eg,
        )

    def condense(
        self,
        target: "EnergyGrid",
        spectrum: np.ndarray,
        /,
        within_group: "WithinGroupSpectrum | None" = None,
        *,
        adjoint_spectrum: np.ndarray | None = None,
    ) -> "Mixture":
        r"""Collapse onto a coarse group structure, spectrum-weighted (rank-0).

        The energy-axis transpose of spatial homogenization, and **data-native** (no
        transport dependency — :class:`~orpheus.data.energy_grid.EnergyGrid`, the overlap
        factory, and the Petrov-Galerkin frame are all data/numerics). Every cross-section
        channel is collapsed from this mixture's fine groups (:attr:`energy_grid`) onto the
        coarser ``target`` structure, preserving each coarse group's reaction rate. The
        fine→coarse fractional-overlap trial is :meth:`EnergyGrid.overlap_to`; the flux
        ``spectrum`` is the TEST weight; the energy measure is counting
        (:meth:`EnergyGrid.as_measure`).

        The collapse has TWO morphisms (the marginalize-vs-average distinction):

        * **average** :math:`G^{-1}M` (``frame.project``) — preserves a reaction RATE:
          the vectors (``SigT``/``SigC``/``SigL``/``SigF``/``SigP``), :math:`\Sigma_G =
          (\sum_g T_{gG}\varphi_g\Sigma_g)/(\sum_g T_{gG}\varphi_g)`; and the **source**
          ``g_from`` axis of each matrix.
        * **marginalize** :math:`@\,T` (the bare overlap table, NO normalisation) —
          preserves MASS: the **sink** ``g_to`` axis of each matrix (``Σ @ T``), and
          :math:`\chi_G = \chi @ T` (a probability over birth groups; preserves
          :math:`\sum\chi=1`). Energy-χ is NOT the weight-1 degenerate of ``project``
          (that would divide by the group count and break the simplex) — it is the
          analysis face *without* the projection normalisation.

        Parameters
        ----------
        target : EnergyGrid
            The coarse target group structure (descending; no more groups than this
            mixture's grid — :meth:`EnergyGrid.overlap_to` enforces the downsampling
            guard).
        spectrum : (ng,) array
            The per-fine-group representative flux :math:`\varphi_g` (the test weight) —
            the solved scalar flux over the cells where the material appears (see
            :meth:`orpheus.sn.solution.Solution.condense`).
        within_group : WithinGroupSpectrum, optional
            The sub-fine-group flux model for straddle apportionment (default 1/E,
            :class:`~orpheus.data.energy_grid.InverseEnergySpectrum`).
        adjoint_spectrum : (ng,) array, optional
            The per-fine-group representative importance :math:`\varphi^*_g` (P6,
            #281).  ``None`` (default) keeps the flux-weighted collapse above,
            bit-identical.  Given, the collapse becomes the **eigenvalue-
            consistent bilinear condensation in the Bell & Glasstone Ch. 6
            convention** (the algebra of record, theorem T6 of
            :mod:`orpheus.derivations.common.homogenization`; memo
            ``p6_literature_memo.md`` Source E): plain condensed-flux carrier,
            flux-weighted-average adjoint carrier :math:`\Psi^{\dagger}_G =
            \langle\varphi^*\varphi\rangle_G/\Phi_G`, and

            * vectors — the diagonal bilinear (B&G (6.135)):
              :math:`\Sigma_{C,G} = \langle\varphi^*\Sigma\varphi\rangle_G/
              \langle\varphi^*\varphi\rangle_G` (the ``average`` morphism with
              the PAIR weight :math:`\varphi^*\!\odot\varphi`);
            * matrices — per-block sink×source (B&G (6.136)): the sink axis
              folds :math:`\varphi^*` and divides by :math:`\Psi^{\dagger}_G`,
              then the source axis flux-averages — the previously weight-free
              ``marginalize`` morphism becomes the adjoint-carried sink
              average;
            * fission — FACTORED (B&G (4.38)+(6.136)): :math:`\nu\Sigma_f`
              flux-weighted, :math:`\chi^{\dagger}_G =
              \langle\chi\varphi^*\rangle_G/\Psi^{\dagger}_G`, with the rank-1
              rescale :math:`s=\sum_G\chi^{\dagger}_G` moved into
              :math:`\nu\Sigma_f` so :math:`\chi_C` stays a simplex
              (:math:`k`-invariant).

            Condensed with the TRUE spectra the bilinear collapse reproduces
            the fine :math:`k` exactly (T6a) and its spectrum-error is second
            order (T6b, B&G (6.90)).

        Returns
        -------
        Mixture
            The condensed (coarse-group) mixture. With ``adjoint_spectrum=None``
            it passes :meth:`assert_balanced` when this one does (every removal
            channel collapses with the same per-coarse-group flux weight
            :math:`\Phi_G`).

            .. warning:: The **bilinear** (``adjoint_spectrum=``) collapse
               does NOT preserve the total-XS balance identity — the classical
               reactivity-vs-rates property of bilinear-weighted constants
               (theorem T4; the sink folds differ per channel).  Do not
               ``assert_balanced`` on a bilinear-condensed mixture.

        Raises
        ------
        ValueError
            If this mixture has no energy grid (``eg is None``), or the ``spectrum``
            length does not match :attr:`ng` (the target/grid group-count mismatch is
            caught by :meth:`EnergyGrid.overlap_to`).
        """
        from orpheus.numerics.basis import WeightedIndicatorBasis
        from orpheus.numerics.frame import PetrovGalerkinFrame

        phi = np.asarray(spectrum, dtype=float)
        if phi.shape != (self.ng,):
            raise ValueError(
                f"spectrum must have shape ({self.ng},), got {phi.shape}."
            )

        source = self.energy_grid                          # raises if eg is None
        trial = source.overlap_to(target, within_group)    # OverlapBasis (the mismatch table)
        table = trial.overlap_table
        frame = PetrovGalerkinFrame(
            trial, source.as_measure(), WeightedIndicatorBasis(trial, phi),
        )

        def average(sig: np.ndarray) -> np.ndarray:
            """Rate-preserving flux-average G⁻¹M over the source-group axis."""
            return frame.project(np.asarray(sig, dtype=float))

        if adjoint_spectrum is None:
            def collapse_matrix(mat: csr_matrix) -> np.ndarray:
                # sink g_to MARGINALIZED (@ T, mass), source g_from AVERAGED (project, rate).
                sink_summed = np.asarray(mat.todense(), dtype=float) @ table
                return average(sink_summed)

            return Mixture.from_dense_channels(
                # vectors → AVERAGE (frame.project = G⁻¹M): rate-preserving flux-average.
                SigC=average(self.SigC),
                SigL=average(self.SigL),
                SigF=average(self.SigF),
                SigP=average(self.SigP),
                SigT=average(self.SigT),
                SigS=[collapse_matrix(s) for s in self.SigS],
                Sig2=collapse_matrix(self.Sig2),
                # χ → MARGINALIZE (bare @ table): a probability over birth groups, mass-
                # preserving (Σχ=1 via the partition-of-unity rows), NOT a flux-average.
                chi=np.asarray(self.chi) @ table,
                eg=target.edges,
            )

        # ── The bilinear (eigenvalue-consistent) collapse — B&G Ch. 6 ──
        phi_star = np.asarray(adjoint_spectrum, dtype=float)
        if phi_star.shape != (self.ng,):
            raise ValueError(
                f"adjoint_spectrum must have shape ({self.ng},), got {phi_star.shape}."
            )

        # The pair frame realizes the diagonal bilinear (B&G 6.135): the
        # ``average`` morphism with the PAIR test weight φ*⊙φ.
        pair_frame = PetrovGalerkinFrame(
            trial, source.as_measure(), WeightedIndicatorBasis(trial, phi_star * phi),
        )

        def bilinear(sig: np.ndarray) -> np.ndarray:
            return pair_frame.project(np.asarray(sig, dtype=float))

        # The adjoint carrier Ψ†_G = ⟨φ*φ⟩_G/Φ_G (B&G 6.126-6.128), with the
        # frames' own Moore-Penrose convention on empty coarse groups.
        pair_G = (phi_star * phi) @ table                  # ⟨φ*φ⟩_G
        flux_G = phi @ table                               # Φ_G
        psi_dag = np.divide(
            pair_G, flux_G, out=np.zeros_like(pair_G), where=flux_G != 0.0,
        )

        def collapse_matrix_bilinear(mat: csr_matrix) -> np.ndarray:
            # sink g_to: the φ*-carried average (fold φ*, divide by Ψ†_G) —
            # the previously weight-free marginalize morphism gains the
            # adjoint carrier (B&G 6.136); source g_from: flux-AVERAGED as
            # always.
            sink_folded = (np.asarray(mat.todense(), dtype=float) * phi_star) @ table
            sink_avg = np.divide(
                sink_folded, psi_dag,
                out=np.zeros_like(sink_folded), where=psi_dag != 0.0,
            )
            return average(sink_avg)

        # The factored fission pair (B&G 4.38 + 6.136): νΣf flux-weighted,
        # χ† adjoint-contracted over Ψ†; the rank-1 rescale s keeps χ a
        # simplex (k-invariant — the dyad's factorization freedom).
        chi_dag = np.divide(
            (np.asarray(self.chi) * phi_star) @ table, psi_dag,
            out=np.zeros_like(psi_dag), where=psi_dag != 0.0,
        )
        s = float(chi_dag.sum())
        chi_c = chi_dag / s if s > 0.0 else chi_dag
        sig_p = average(self.SigP) * (s if s > 0.0 else 1.0)

        return Mixture.from_dense_channels(
            SigC=bilinear(self.SigC),
            SigL=bilinear(self.SigL),
            SigF=bilinear(self.SigF),
            SigP=sig_p,
            SigT=bilinear(self.SigT),
            SigS=[collapse_matrix_bilinear(m) for m in self.SigS],
            Sig2=collapse_matrix_bilinear(self.Sig2),
            chi=chi_c,
            eg=target.edges,
        )


def production_weighted_chi(
    isotopes: list[Isotope],
    sigF: np.ndarray,
    aDen: np.ndarray,
    fissile_indices: list[int],
) -> np.ndarray:
    r"""Mixture fission spectrum: the production-weighted convex average of
    the fissile isotopes' emission spectra.

    .. math:: \chi_{\rm mix} = \sum_i w_i\,\chi_i,\quad
              w_i = \frac{a_i \sum_g \bar\nu_{i,g}\sigma_{f,i,g}}
                         {\sum_j a_j \sum_g \bar\nu_{j,g}\sigma_{f,j,g}}

    where the per-isotope production rate :math:`p_i = a_i \sum_g
    \bar\nu_{i,g}\sigma_{f,i,g}` is the flat-flux representative fission
    production density (the standard data-prep weighting, :math:`\phi_g
    \equiv 1`).

    HONEST SCOPE (vv Mode-7). The truly-correct mixture spectrum is
    FLUX-weighted (:math:`w_i \propto \sum_g \nu\Sigma_{f,i,g}\,\phi_g`); a
    static data-prep value cannot capture :math:`\phi`. This is a
    DOCUMENTED approximation, NOT a separate issue — the flat-flux
    production weighting is the conventional data-prep choice. A convex
    combination of simplices is a simplex, so the result passes
    :meth:`~orpheus.data.emission_spectrum.EmissionSpectrum.assert_normalized`
    when the per-isotope :math:`\chi_i` are themselves simplices.

    Parameters
    ----------
    isotopes : list[Isotope]
        Microscopic cross-section data; ``isotopes[i].nubar`` (NG,) and
        ``isotopes[i].chi`` (NG, simplex) are read.
    sigF : (n_iso, NG) array
        Per-isotope fission XS interpolated at the converged sigma-zeros
        (the same array :func:`compute_macro_xs` builds).
    aDen : (n_iso,) array
        Number densities in 1/(barn*cm).
    fissile_indices : list[int]
        Indices of the fissile isotopes to average over (MUST be non-empty;
        the caller guards the empty case).

    Returns
    -------
    (NG,) array
        The production-weighted convex average emission spectrum. Returned
        as a plain ``np.ndarray``; the :class:`Mixture` constructor coerces
        it to a validated :class:`EmissionSpectrum` in ``__post_init__``.
    """
    production = np.array(
        [aDen[i] * float(np.sum(isotopes[i].nubar * sigF[i])) for i in fissile_indices]
    )
    total = production.sum()
    # fissile_indices is selected by is_fissile (σf>0) and real isotopes
    # carry ν̄>0 where σf>0, so total>0 in practice. Guard the degenerate
    # ν̄≡0 (zero total production) case rather than dividing by zero.
    if total > 0:
        weights = production / total
    else:
        weights = np.full(len(fissile_indices), 1.0 / len(fissile_indices))
    # The convex average χ_mix = Σ_i w_i χ_i as a weights·spectra contraction:
    # stack the per-isotope simplices into (n_fissile, NG), then weights (1, n)
    # @ spectra (n, NG) → (NG,). Reads like the math and types as an ndarray.
    fissile_spectra = np.array([np.asarray(isotopes[i].chi) for i in fissile_indices])
    return weights @ fissile_spectra


def compute_macro_xs(
    isotopes: list[Isotope],
    number_densities: np.ndarray,
    escape_xs: float = 0.0,
    n_legendre: int = 3,
    fissile_indices: Optional[list[int]] = None,
) -> Mixture:
    """Compute macroscopic cross sections for a mixture of isotopes.

    Parameters
    ----------
    isotopes : list[Isotope]
        Microscopic cross section data for each isotope.
    number_densities : (n_iso,) array
        Number densities in 1/(barn*cm).
    escape_xs : float
        Escape cross section in 1/cm (0 for infinite medium).
    n_legendre : int
        Number of Legendre scattering components (default 3).
    fissile_indices : list[int] or None
        Indices into `isotopes` for fissile nuclides.  If None, auto-detected.

    Returns
    -------
    Mixture with all macroscopic cross sections.
    """
    n_iso = len(isotopes)
    aDen = np.asarray(number_densities)
    eg = isotopes[0].eg

    # --- Sigma-zero iterations ---
    print("  Sigma-zero iterations...", end=" ", flush=True)
    sig0 = solve_sigma_zeros(isotopes, aDen, escape_xs)
    print("done.")

    # --- Interpolate microscopic XS at converged sigma-zeros ---
    print("  Interpolating cross sections...", end=" ", flush=True)

    sigC = np.array([interp_xs_field(iso.sigC, iso, sig0[i]) for i, iso in enumerate(isotopes)])
    sigL = np.array([interp_xs_field(iso.sigL, iso, sig0[i]) for i, iso in enumerate(isotopes)])
    sigF = np.array([interp_xs_field(iso.sigF, iso, sig0[i]) for i, iso in enumerate(isotopes)])

    sigS_list: list[list[csr_matrix]] = []
    for j in range(n_legendre):
        sigS_j = [interp_sig_s(iso, j, sig0[i]) for i, iso in enumerate(isotopes)]
        sigS_list.append(sigS_j)

    print("done.")

    # --- Sum to macroscopic XS ---
    # Vector XS: (NG,) = sum_i micro_i(NG) * N_i
    SigC = sigC.T @ aDen
    SigL = sigL.T @ aDen
    SigF = sigF.T @ aDen

    # Production XS: only from fissile isotopes. The start value is the
    # empty-sum identity (zero production for a non-fissile mixture).
    if fissile_indices is None:
        fissile_indices = [i for i, iso in enumerate(isotopes) if iso.is_fissile]
    SigP = sum(
        (isotopes[i].nubar * sigF[i] * aDen[i] for i in fissile_indices),
        start=np.zeros(NG),
    )

    # Scattering matrices
    SigS = [
        sum(
            (sigS_list[j][i] * aDen[i] for i in range(n_iso)),
            start=csr_matrix((NG, NG)),
        )
        for j in range(n_legendre)
    ]

    # (n,2n) matrix
    Sig2 = sum(
        (iso.sig2 * aDen[i] for i, iso in enumerate(isotopes)),
        start=csr_matrix((NG, NG)),
    )

    # Total XS
    SigT = SigC + SigL + SigF + np.array(SigS[0].sum(axis=1)).ravel() + np.array(Sig2.sum(axis=1)).ravel()

    # Fission spectrum — production-weighted convex average over ALL fissile
    # isotopes (flat-flux representative weighting; see
    # `production_weighted_chi`). The Mixture constructor coerces this plain
    # array to a validated EmissionSpectrum.
    chi = np.zeros(NG)
    if fissile_indices:
        chi = production_weighted_chi(isotopes, sigF, aDen, fissile_indices)

    mix = Mixture(
        SigC=SigC, SigL=SigL, SigF=SigF, SigP=SigP, SigT=SigT,
        SigS=SigS, Sig2=Sig2, chi=chi, eg=eg,
    )
    # Real-path guard: SigT was just DERIVED via the balance identity, so a
    # physical mixture always balances. This is a free regression guard
    # against a future derivation bug (it catches, it does not break).
    mix.assert_balanced()
    return mix
