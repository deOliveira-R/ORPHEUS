"""Macroscopic cross section computation for isotope mixtures.

A Mixture holds the isotopes, their number densities, and the resulting
macroscopic cross sections used by reactor solvers.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from scipy.sparse import csr_matrix

from orpheus.data.emission_spectrum import enforce_emission_spectrum
from orpheus.data.micro_xs.isotope import NG, Isotope
from .interpolation import interp_sig_s, interp_xs_field
from .sigma_zeros import solve_sigma_zeros


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
        (Galerkin projection from fine to coarse groups). That step
        depends on the grid; it cannot run on synthetic mixtures.
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

    # Production XS: only from fissile isotopes
    if fissile_indices is None:
        fissile_indices = [i for i, iso in enumerate(isotopes) if iso.is_fissile]
    SigP = sum(
        isotopes[i].nubar * sigF[i] * aDen[i] for i in fissile_indices
    ) if fissile_indices else np.zeros(NG)

    # Scattering matrices
    SigS = []
    for j in range(n_legendre):
        mat = sum(sigS_list[j][i] * aDen[i] for i in range(n_iso))
        SigS.append(mat)

    # (n,2n) matrix
    Sig2 = sum(iso.sig2 * aDen[i] for i, iso in enumerate(isotopes))

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
