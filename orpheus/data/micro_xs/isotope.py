"""Data model for microscopic cross section data of a single isotope."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.sparse import csr_matrix

from orpheus.data.emission_spectrum import enforce_emission_spectrum

NG = 421  # number of energy groups


@dataclass
class Isotope:
    """Microscopic cross section data for one isotope at one temperature.

    All cross sections are in barns; energy boundaries in eV.
    Dimensions:
        sig0    : (n_sig0,)           — background cross section base points
        sigC    : (n_sig0, NG)        — radiative capture
        sigL    : (n_sig0, NG)        — (n,alpha)
        sigF    : (n_sig0, NG)        — fission
        sigT    : (n_sig0, NG)        — total
        sigS    : [3][n_sig0] of (NG, NG) sparse  — scattering (3 Legendre)
        sig2    : (NG, NG) sparse     — (n,2n) matrix
        nubar   : (NG,)              — average neutrons per fission
        chi     : (NG,)              — fission spectrum
        eg      : (NG+1,)            — energy group boundaries
    """

    name: str
    aw: float  # atomic weight (amu)
    temp: float  # temperature (K)
    eg: np.ndarray  # energy group boundaries
    sig0: np.ndarray  # sigma-zero base points

    sigC: np.ndarray
    sigL: np.ndarray
    sigF: np.ndarray
    sigT: np.ndarray

    nubar: np.ndarray  # (NG,)
    chi: np.ndarray  # (NG,)

    sigS: list[list[csr_matrix]] = field(default_factory=list)  # [legendre][sig0_idx]
    sig2: csr_matrix = field(default_factory=lambda: csr_matrix((NG, NG)))

    def __post_init__(self) -> None:
        # Coerce chi to the validated value-object and enforce the simplex
        # / null law at the data source. χ is consumed only as a fission
        # SOURCE (χ·νΣ_f·φ), so the law keys on PRODUCTION (νΣ_f > 0): a
        # producing isotope's spectrum is a probability simplex; a
        # non-producing isotope emits no fission neutrons, so its spectrum
        # is identically zero.
        self.chi = enforce_emission_spectrum(self.chi, is_producing=self.is_producing)

    @property
    def n_sig0(self) -> int:
        return len(self.sig0)

    @property
    def ng(self) -> int:
        return len(self.eg) - 1

    @property
    def is_fissile(self) -> bool:
        """``True`` iff the isotope can fission (nonzero fission XS in any group).

        Distinct from :attr:`is_producing`: ``is_fissile`` is the
        cross-section question ("can it fission?") consumed by
        ``compute_macro_xs`` (``fissile_indices``); ``is_producing`` is the
        emission question ("does it emit fission neutrons?", ``νΣ_f > 0``)
        that the χ simplex/null law keys on.
        """
        return np.any(self.sigF > 0)

    @property
    def is_producing(self) -> bool:
        r"""``True`` iff the isotope emits fission neutrons (:math:`\nu\Sigma_f > 0`).

        Broadcast over the ``(n_sig0,)`` background-XS base points:
        ``nubar`` is ``(NG,)`` and ``sigF`` is ``(n_sig0, NG)``, so
        ``nubar * sigF`` is ``(n_sig0, NG)``. This is the predicate the χ
        emission-spectrum law keys on (a valid simplex is required exactly
        where production is nonzero).
        """
        return bool(np.any(self.nubar * self.sigF > 0))
