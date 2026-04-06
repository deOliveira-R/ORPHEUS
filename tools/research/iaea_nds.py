"""IAEA Nuclear Data Services (NDS) LiveChart API client.

Provides typed access to the IAEA NDS LiveChart of Nuclides REST API
for querying nuclear structure data: ground state properties, gamma
transitions, and energy levels.

This is nuclear DATA, not literature — use it for physical constants,
half-lives, binding energies, decay modes, and gamma energies.

No authentication required. Returns CSV parsed into dataclasses.

Available data via this API:

- ``ground_states`` — mass, binding energy, half-life, spin-parity, abundance
- ``levels`` — excited nuclear energy levels
- ``gammas`` — gamma transition energies, intensities, multipolarities

Note: fission yields are NOT available through this endpoint.
Use ENDF libraries for fission yield data.

Usage::

    from tools.research.iaea_nds import (
        ground_state, gammas, levels,
    )

    # Ground state properties of U-235
    gs = ground_state("U235")
    print(gs.half_life, gs.binding_energy, gs.spin_parity)

    # Gamma transitions of Co-60
    gs_gammas = gammas("Co60")
    for g in gs_gammas[:5]:
        print(f"{g.energy} keV, I={g.relative_intensity}%, {g.multipolarity}")

    # All uranium isotopes
    from tools.research.iaea_nds import ground_states_by_z
    isotopes = ground_states_by_z(92)
    for iso in isotopes:
        print(iso.summary())
"""

from __future__ import annotations

import csv
import io
import time
from dataclasses import dataclass
from typing import Optional

import requests

BASE_URL = "https://nds.iaea.org/relnsd/v1/data"

_MIN_REQUEST_INTERVAL = 0.5
_last_request_time = 0.0


def _throttle() -> None:
    global _last_request_time
    elapsed = time.monotonic() - _last_request_time
    if elapsed < _MIN_REQUEST_INTERVAL:
        time.sleep(_MIN_REQUEST_INTERVAL - elapsed)
    _last_request_time = time.monotonic()


def _fetch_csv(params: dict) -> list[dict]:
    """Fetch CSV data from the NDS API and parse into dicts."""
    _throttle()
    resp = requests.get(BASE_URL, params=params, timeout=30)
    resp.raise_for_status()
    text = resp.text.strip()
    if not text or "," not in text:
        # API returns a bare number (error code) on invalid field names.
        return []
    reader = csv.DictReader(io.StringIO(text))
    return list(reader)


def _parse_float(value: str) -> Optional[float]:
    """Parse a float from CSV, returning None for empty/invalid."""
    if not value or value.strip() == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


# ── Data classes ──────────────────────────────────────────────────────


@dataclass
class GroundState:
    """Ground state nuclear properties of a nuclide.

    CSV columns: z, n, symbol, radius, unc_r, abundance, unc_a,
    energy_shift, energy, unc_e, ripl_shift, jp, half_life, operator_hl,
    unc_hl, unit_hl, half_life_sec, unc_hls, decay_1, decay_1_%,
    unc_d1, decay_2, decay_2_%, unc_d2, decay_3, decay_3_%, unc_d3,
    binding, unc_b, atomic_mass, unc_am, mass_excess, unc_me, ...
    """

    z: int
    n: int
    symbol: str
    mass_number: int
    atomic_mass: Optional[float]  # u
    mass_excess: Optional[float]  # keV
    binding_energy: Optional[float]  # keV/nucleon
    half_life: str
    half_life_seconds: Optional[float]
    decay_modes: str
    spin_parity: str
    abundance: Optional[float]  # natural isotopic abundance (%)
    radius: Optional[float]  # fm

    @classmethod
    def from_csv(cls, row: dict) -> GroundState:
        z = int(row.get("z", 0))
        n = int(row.get("n", 0))
        return cls(
            z=z,
            n=n,
            symbol=row.get("symbol", ""),
            mass_number=z + n,
            atomic_mass=_parse_float(row.get("atomic_mass", "")),
            mass_excess=_parse_float(row.get("mass_excess", "")),
            binding_energy=_parse_float(row.get("binding", "")),
            half_life=row.get("half_life", ""),
            half_life_seconds=_parse_float(row.get("half_life_sec", "")),
            decay_modes=row.get("decay_1", ""),
            spin_parity=row.get("jp", ""),
            abundance=_parse_float(row.get("abundance", "")),
            radius=_parse_float(row.get("radius", "")),
        )

    def summary(self) -> str:
        return (
            f"{self.symbol}-{self.mass_number} "
            f"(Z={self.z}, N={self.n}): "
            f"t½={self.half_life}, "
            f"BE={self.binding_energy} keV/A, "
            f"Jπ={self.spin_parity}"
        )


@dataclass
class GammaTransition:
    """A gamma transition between nuclear levels (from ``fields=gammas``).

    CSV columns: z, n, symbol, start_level_idx, start_level_energy,
    unc_sle, start_level_jp, end_level_idx, end_level_energy, unc_ele,
    end_level_jp, gamma_idx, energy, unc_en, relative_intensity, unc_ri,
    multipolarity, mixing_ratio, unc_mr, ...
    """

    z: int
    n: int
    symbol: str
    energy: Optional[float]  # keV
    energy_uncertainty: Optional[float]
    relative_intensity: Optional[float]  # %
    intensity_uncertainty: Optional[float]
    multipolarity: str
    start_level_energy: Optional[float]  # keV
    start_level_jp: str
    end_level_energy: Optional[float]  # keV
    end_level_jp: str

    @classmethod
    def from_csv(cls, row: dict) -> GammaTransition:
        return cls(
            z=int(row.get("z", 0)),
            n=int(row.get("n", 0)),
            symbol=row.get("symbol", ""),
            energy=_parse_float(row.get("energy", "")),
            energy_uncertainty=_parse_float(row.get("unc_en", "")),
            relative_intensity=_parse_float(row.get("relative_intensity", "")),
            intensity_uncertainty=_parse_float(row.get("unc_ri", "")),
            multipolarity=row.get("multipolarity", ""),
            start_level_energy=_parse_float(row.get("start_level_energy", "")),
            start_level_jp=row.get("start_level_jp", ""),
            end_level_energy=_parse_float(row.get("end_level_energy", "")),
            end_level_jp=row.get("end_level_jp", ""),
        )

    def summary(self) -> str:
        e = f"{self.energy:.2f}" if self.energy is not None else "?"
        i = f"{self.relative_intensity:.1f}" if self.relative_intensity is not None else "?"
        return f"{e} keV ({self.start_level_jp} → {self.end_level_jp}), I={i}%, {self.multipolarity}"


@dataclass
class DecayRadiation:
    """A decay radiation entry (from ``fields=decay_rads``).

    CSV columns: energy, unc_en, intensity, unc_i, start_level_hl,
    start_level_energy, end_level_hl, end_level_energy, multipolarity,
    mixing_ratio, unc_mr, conversion_coeff, unc_cc, p_z, p_n, p_symbol, ...
    """

    parent_z: int
    parent_n: int
    parent_symbol: str
    energy: Optional[float]  # keV
    energy_uncertainty: Optional[float]
    intensity: Optional[float]  # per 100 decays
    intensity_uncertainty: Optional[float]
    multipolarity: str
    start_level_energy: Optional[float]
    end_level_energy: Optional[float]

    @classmethod
    def from_csv(cls, row: dict) -> DecayRadiation:
        return cls(
            parent_z=int(row.get("p_z", 0)),
            parent_n=int(row.get("p_n", 0)),
            parent_symbol=row.get("p_symbol", ""),
            energy=_parse_float(row.get("energy", "")),
            energy_uncertainty=_parse_float(row.get("unc_en", "")),
            intensity=_parse_float(row.get("intensity", "")),
            intensity_uncertainty=_parse_float(row.get("unc_i", "")),
            multipolarity=row.get("multipolarity", ""),
            start_level_energy=_parse_float(row.get("start_level_energy", "")),
            end_level_energy=_parse_float(row.get("end_level_energy", "")),
        )

    def summary(self) -> str:
        e = f"{self.energy:.2f}" if self.energy is not None else "?"
        i = f"{self.intensity:.2f}" if self.intensity is not None else "?"
        return f"{e} keV, I={i}%, {self.multipolarity}"


@dataclass
class NuclearLevel:
    """An excited nuclear level.

    CSV columns: z, n, symbol, idx, energy_shift, energy, unc_e,
    ripl_shift, jp, jp_order, half_life, operator_hl, unc_hl, unit_hl,
    half_life_sec, unc_hls, ...
    """

    z: int
    n: int
    symbol: str
    energy: Optional[float]  # keV
    energy_uncertainty: Optional[float]
    spin_parity: str
    half_life: str
    half_life_seconds: Optional[float]

    @classmethod
    def from_csv(cls, row: dict) -> NuclearLevel:
        return cls(
            z=int(row.get("z", 0)),
            n=int(row.get("n", 0)),
            symbol=row.get("symbol", ""),
            energy=_parse_float(row.get("energy", "")),
            energy_uncertainty=_parse_float(row.get("unc_e", "")),
            spin_parity=row.get("jp", ""),
            half_life=row.get("half_life", ""),
            half_life_seconds=_parse_float(row.get("half_life_sec", "")),
        )

    def summary(self) -> str:
        e = f"{self.energy:.2f}" if self.energy is not None else "?"
        return f"E={e} keV, Jπ={self.spin_parity}, t½={self.half_life}"


# ── API functions ─────────────────────────────────────────────────────


def ground_state(nuclide: str) -> GroundState:
    """Get ground state properties of a nuclide.

    Parameters
    ----------
    nuclide : str
        Nuclide identifier (e.g. "U235", "Pu239", "Co60", "H1").

    Returns
    -------
    GroundState
    """
    rows = _fetch_csv({"fields": "ground_states", "nuclides": nuclide})
    if not rows:
        raise ValueError(f"No ground state data for {nuclide}")
    return GroundState.from_csv(rows[0])


def ground_states(nuclides: list[str]) -> list[GroundState]:
    """Get ground state properties for multiple nuclides.

    Parameters
    ----------
    nuclides : list[str]
        Nuclide identifiers (e.g. ["U235", "U238", "Pu239"]).

    Returns
    -------
    list[GroundState]
    """
    rows = _fetch_csv({
        "fields": "ground_states",
        "nuclides": ",".join(nuclides),
    })
    return [GroundState.from_csv(row) for row in rows]


def ground_states_by_z(z: int) -> list[GroundState]:
    """Get ground states for all isotopes of an element.

    Parameters
    ----------
    z : int
        Atomic number.

    Returns
    -------
    list[GroundState]
    """
    rows = _fetch_csv({"fields": "ground_states", "nuclides": f"z={z}"})
    return [GroundState.from_csv(row) for row in rows]


def gammas(nuclide: str) -> list[GammaTransition]:
    """Get gamma transitions for a nuclide.

    Returns gamma-ray energies, relative intensities, and multipolarities
    from excited nuclear levels.

    Parameters
    ----------
    nuclide : str
        Nuclide identifier (e.g. "Co60", "Cs137", "Fe56").

    Returns
    -------
    list[GammaTransition]
    """
    rows = _fetch_csv({"fields": "gammas", "nuclides": nuclide})
    return [GammaTransition.from_csv(row) for row in rows]


def decay_radiations(
    nuclide: str,
    rad_types: str = "g",
) -> list[DecayRadiation]:
    """Get decay radiation data for a nuclide.

    Parameters
    ----------
    nuclide : str
        Nuclide identifier (e.g. "Co60", "Cs137").
    rad_types : str
        Radiation types to include (comma-separated):
        "g" (gamma), "b-" (beta-minus), "bp" (beta-plus),
        "a" (alpha), "x" (X-ray), "e" (conversion electron).
        Default "g" (gamma only). Use "g,b-,a" for multiple.

    Returns
    -------
    list[DecayRadiation]
    """
    rows = _fetch_csv({
        "fields": "decay_rads",
        "nuclides": nuclide,
        "rad_types": rad_types,
    })
    return [DecayRadiation.from_csv(row) for row in rows]


def levels(nuclide: str) -> list[NuclearLevel]:
    """Get excited nuclear energy levels for a nuclide.

    Parameters
    ----------
    nuclide : str
        Nuclide identifier (e.g. "Fe56", "U238").

    Returns
    -------
    list[NuclearLevel]
    """
    rows = _fetch_csv({"fields": "levels", "nuclides": nuclide})
    return [NuclearLevel.from_csv(row) for row in rows]
