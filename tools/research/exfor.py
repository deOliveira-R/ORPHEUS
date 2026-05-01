"""EXFOR (Experimental Nuclear Reaction Database) client.

Programmatic access to the IAEA Nuclear Data Section's experimental
reaction database. EXFOR contains compiled measurements of cross
sections, angular distributions, secondary particle spectra, and
fission yields from nuclear physics experiments going back to the
1930s. It is the canonical source of measured data underlying every
modern evaluated nuclear data library (ENDF/B, JEFF, JENDL, BROND,
CENDL).

This client uses the IAEA-NDS open-data GitHub repositories rather
than the (poorly documented) Cloudflare-fronted Dash REST API:

- ``IAEA-NDS/exfortables_py`` — pre-tabulated ``(x, dx, y, dy)``
  text files organised as ``{projectile}/{target}/{reaction}/
  {observable}/{filename}.txt``. One file per EXFOR subentry. This
  is the next-generation rebuild of the classic EXFORTABLES (TALYS).
- ``IAEA-NDS/exfor_json`` — full ENTRY/SUBENTRY JSON conversions of
  the EXFOR master files; useful when you need the complete
  bibliographic record, experimental conditions, or references.

GitHub raw content has effectively no rate limit; the GitHub
contents API (used for directory listing) is throttled at 60/h
unauthenticated, 5000/h with a token. Set ``GITHUB_TOKEN`` to lift.

Usage::

    from tools.research.exfor import (
        list_datasets, get_dataset, get_entry_json, parse_filename,
    )

    # 1. Discover available datasets for a reaction
    files = list_datasets(projectile="n", target="U-235",
                          reaction="n-f", observable="xs")
    # [DatasetRef(filename='U-235_n-f_A.D.Carlson-14015-002-0-1991.txt', ...)]

    # 2. Pull the actual data points
    ds = get_dataset(files[0])
    ds.x, ds.y, ds.dx, ds.dy        # numpy-friendly lists of floats
    ds.author, ds.year, ds.entry_id, ds.reaction
    ds.facility, ds.institute, ds.reference

    # 3. Get the full JSON entry record (BIB block + every subentry)
    entry = get_entry_json("14015")
    entry["bib_record"]["AUTHOR"]
    entry["data_tables"]            # all subentries
"""

from __future__ import annotations

import os
import re
import time
from dataclasses import dataclass, field
from typing import Optional, Union

import requests

# GitHub-hosted EXFOR data repositories (IAEA-NDS, MIT licence).
_EXFORTABLES_OWNER = "IAEA-NDS"
_EXFORTABLES_REPO = "exfortables_py"
_EXFOR_JSON_OWNER = "IAEA-NDS"
_EXFOR_JSON_REPO = "exfor_json"
_GH_API = "https://api.github.com"
_GH_RAW = "https://raw.githubusercontent.com"

GITHUB_TOKEN: Optional[str] = os.environ.get("GITHUB_TOKEN")

# Reaction-code shorthands for common neutron-induced reactions.
# Maps a friendly form to the directory name used in exfortables_py.
COMMON_REACTIONS = {
    # Cross sections
    "(n,f)": "n-f",
    "(n,fission)": "n-f",
    "(n,2n)": "n-2n",
    "(n,3n)": "n-3n",
    "(n,4n)": "n-4n",
    "(n,gamma)": "n-g",
    "(n,g)": "n-g",
    "(n,total)": "n-tot",
    "(n,tot)": "n-tot",
    "(n,el)": "n-el",
    "(n,elastic)": "n-el",
    "(n,inl)": "n-inl",
    "(n,inelastic)": "n-inl",
    "(n,non)": "n-non",
    "(n,abs)": "n-abs",
    "(n,p)": "n-p",
    "(n,a)": "n-a",
    "(n,alpha)": "n-a",
    "(n,scat)": "n-sct",
    "(n,x)": "n-x",
}

# Reaction-MT mapping (subset relevant to reactor physics). MT codes
# follow ENDF-6 conventions.
MT_MAP = {
    "n-tot": 1,
    "n-el": 2,
    "n-non": 3,
    "n-inl": 4,
    "n-2n": 16,
    "n-3n": 17,
    "n-f": 18,
    "n-g": 102,
    "n-p": 103,
    "n-a": 107,
}

_MIN_REQUEST_INTERVAL = 1.0  # polite throttle for GitHub
_last_request_time = 0.0


def _throttle() -> None:
    global _last_request_time
    elapsed = time.monotonic() - _last_request_time
    if elapsed < _MIN_REQUEST_INTERVAL:
        time.sleep(_MIN_REQUEST_INTERVAL - elapsed)
    _last_request_time = time.monotonic()


def _gh_headers() -> dict:
    h = {"Accept": "application/vnd.github+json"}
    if GITHUB_TOKEN:
        h["Authorization"] = f"Bearer {GITHUB_TOKEN}"
    return h


def normalize_reaction(reaction: str) -> str:
    """Convert a friendly reaction code into the exfortables_py directory form.

    Accepts ``"(n,f)"``, ``"n,f"``, ``"n-f"``, ``"N,F"`` etc.
    """
    r = reaction.strip().lower()
    if r in COMMON_REACTIONS:
        return COMMON_REACTIONS[r]
    # Strip parens, normalise separators.
    r = r.strip("()").replace(",", "-").replace(" ", "")
    return r


def normalize_target(target: str) -> str:
    """Normalise a target nuclide string to the ``El-A`` form.

    Accepts ``"U235"``, ``"u-235"``, ``"U-235"``, ``"235U"``,
    ``"235-U"``. Returns ``"U-235"`` (capitalised symbol, dash, mass).
    """
    t = target.strip()
    m = re.match(r"^\s*(\d+)\s*-?\s*([A-Za-z]{1,3})\s*$", t)  # 235U / 235-U
    if m:
        a, sym = m.groups()
        return f"{sym.capitalize()}-{a}"
    m = re.match(r"^\s*([A-Za-z]{1,3})\s*-?\s*(\d+)\s*$", t)  # U235 / U-235
    if m:
        sym, a = m.groups()
        return f"{sym.capitalize()}-{a}"
    return t  # leave unrecognised input alone (e.g. natural abundance "Ag-0")


# ── Data classes ──────────────────────────────────────────────────────


@dataclass
class DatasetRef:
    """A pointer to a tabulated EXFOR dataset (one subentry)."""

    filename: str
    url: str  # raw GitHub URL — pass directly to get_dataset()
    projectile: str = ""
    target: str = ""
    reaction: str = ""
    observable: str = ""

    # Parsed from filename
    author: str = ""
    entry_id: str = ""    # e.g. "14015"
    subentry: str = ""    # e.g. "002"
    pointer: str = ""     # e.g. "0" or "1"
    year: str = ""

    @classmethod
    def from_github_item(
        cls,
        item: dict,
        *,
        projectile: str,
        target: str,
        reaction: str,
        observable: str,
    ) -> DatasetRef:
        name = item.get("name", "")
        url = item.get("download_url") or ""
        ref = cls(
            filename=name,
            url=url,
            projectile=projectile,
            target=target,
            reaction=reaction,
            observable=observable,
        )
        parsed = parse_filename(name)
        ref.author = parsed.get("author", "")
        ref.entry_id = parsed.get("entry_id", "")
        ref.subentry = parsed.get("subentry", "")
        ref.pointer = parsed.get("pointer", "")
        ref.year = parsed.get("year", "")
        return ref

    def summary(self) -> str:
        return (
            f"[EXFOR:{self.entry_id}-{self.subentry}-{self.pointer}] "
            f"{self.author or 'Unknown'} ({self.year or 'n.d.'}): "
            f"{self.target} {self.reaction} {self.observable}"
        )


@dataclass
class Dataset:
    """A parsed EXFOR tabulated dataset (one subentry) + provenance."""

    # Bibliographic / experimental conditions
    title: str = ""
    target: str = ""
    target_z: str = ""
    target_a: str = ""
    process: str = ""               # e.g. "N,F"
    mf: str = ""                    # ENDF MF
    mt: str = ""                    # ENDF MT
    energy_range: str = ""
    residual: str = ""
    entry_id: str = ""              # "entry-subentry-pointer"
    reaction_code: str = ""         # e.g. "(92-U-235(N,F),,SIG)"
    author: str = ""
    institute: str = ""
    reference: str = ""
    year: str = ""
    facility: str = ""
    master_file_url: str = ""       # link to the .x4 EXFOR master
    nds_url: str = ""               # IAEA NDS landing page

    # Data column labels (parsed from the data header line)
    column_labels: list[str] = field(default_factory=list)

    # Tabular data — same order as column_labels
    rows: list[list[float]] = field(default_factory=list)

    @property
    def x(self) -> list[float]:
        """First column (typically incident energy)."""
        return [r[0] for r in self.rows] if self.rows else []

    @property
    def dx(self) -> list[float]:
        """Second column (typically energy uncertainty)."""
        return [r[1] for r in self.rows] if self.rows and len(self.rows[0]) > 1 else []

    @property
    def y(self) -> list[float]:
        """Third column (typically the measured observable)."""
        return [r[2] for r in self.rows] if self.rows and len(self.rows[0]) > 2 else []

    @property
    def dy(self) -> list[float]:
        """Fourth column (typically observable uncertainty)."""
        return [r[3] for r in self.rows] if self.rows and len(self.rows[0]) > 3 else []

    @property
    def n_points(self) -> int:
        return len(self.rows)

    def summary(self) -> str:
        return (
            f"[EXFOR:{self.entry_id}] {self.author or '?'} "
            f"({self.year or 'n.d.'}): {self.title or self.reaction_code} "
            f"— {self.n_points} points"
        )


# ── Filename parser ───────────────────────────────────────────────────


_FILENAME_RE = re.compile(
    r"^(?P<target>[^_]+)_(?P<reaction>[^_]+)_(?P<author>.+?)-"
    r"(?P<entry>[A-Z]?\d+)-(?P<subent>\d+)-(?P<pointer>\d+)-"
    r"(?P<year>\d{4})\.txt$"
)


def parse_filename(filename: str) -> dict:
    """Parse an exfortables_py filename into its components.

    Filename pattern (one EXFOR subentry per file):
    ``{target}_{reaction}_{author}-{entry_id}-{subentry}-{pointer}-{year}.txt``

    Example
    -------
    >>> parse_filename("U-235_n-f_A.D.Carlson-14015-002-0-1991.txt")
    {'target': 'U-235', 'reaction': 'n-f', 'author': 'A.D.Carlson',
     'entry_id': '14015', 'subentry': '002', 'pointer': '0', 'year': '1991'}
    """
    m = _FILENAME_RE.match(filename)
    if not m:
        return {}
    d = m.groupdict()
    # Rename to the names used by DatasetRef.
    return {
        "target": d["target"],
        "reaction": d["reaction"],
        "author": d["author"],
        "entry_id": d["entry"],
        "subentry": d["subent"],
        "pointer": d["pointer"],
        "year": d["year"],
    }


# ── Discovery ─────────────────────────────────────────────────────────


def list_datasets(
    *,
    projectile: str,
    target: str,
    reaction: str,
    observable: str = "xs",
) -> list[DatasetRef]:
    """List every EXFOR dataset in IAEA-NDS/exfortables_py for a reaction.

    Parameters
    ----------
    projectile : str
        ``"n"`` (neutron), ``"p"`` (proton), ``"g"`` (gamma),
        ``"a"`` (alpha), ``"d"`` (deuteron), ``"t"`` (triton),
        ``"h"`` (helium-3), ``"0"`` (spontaneous reaction with no
        projectile, e.g. spontaneous fission of Cm-244).
    target : str
        Target nuclide. Accepts loose forms — :func:`normalize_target`
        coerces ``"U235"`` and ``"235U"`` into ``"U-235"``. Use
        ``"El-0"`` for natural abundance (``"Ag-0"``).
    reaction : str
        Reaction code. Accepts ``"(n,f)"``, ``"n-f"``, ``"n,f"``.
        :func:`normalize_reaction` recognises the common variants;
        otherwise pass the directory form directly. See
        :data:`COMMON_REACTIONS`.
    observable : str
        Observable type. Defaults to ``"xs"`` (cross section). Other
        common values: ``"angle"`` (angular distribution),
        ``"energy"`` (secondary spectrum), ``"fission/yield/cumulative"``,
        ``"fission/yield/independent"``, ``"fission/yield/primary"``.

    Returns
    -------
    list[DatasetRef]
        One entry per EXFOR subentry. Empty list if the reaction has
        no measurements deposited.

    Notes
    -----
    Hits the GitHub contents API. With no ``GITHUB_TOKEN`` set,
    you have 60 anonymous requests per hour — fine for interactive
    use but not for batch enumeration of many reactions.
    """
    tgt = normalize_target(target)
    rxn = normalize_reaction(reaction)
    path = f"{projectile}/{tgt}/{rxn}/{observable}"

    url = (
        f"{_GH_API}/repos/{_EXFORTABLES_OWNER}/{_EXFORTABLES_REPO}/"
        f"contents/{path}?ref=main"
    )
    _throttle()
    resp = requests.get(url, headers=_gh_headers(), timeout=30)
    if resp.status_code == 404:
        return []
    resp.raise_for_status()
    items = resp.json()
    if not isinstance(items, list):
        return []

    return [
        DatasetRef.from_github_item(
            item,
            projectile=projectile,
            target=tgt,
            reaction=rxn,
            observable=observable,
        )
        for item in items
        if item.get("type") == "file" and item.get("name", "").endswith(".txt")
    ]


# ── Dataset retrieval ─────────────────────────────────────────────────


def get_dataset(ref_or_url: Union[DatasetRef, str]) -> Dataset:
    """Fetch and parse a single EXFOR tabulated dataset.

    Parameters
    ----------
    ref_or_url : DatasetRef or str
        Either a :class:`DatasetRef` returned by :func:`list_datasets`,
        or a raw GitHub URL pointing at a ``.txt`` file in the
        ``exfortables_py`` repo.

    Returns
    -------
    Dataset
    """
    url = ref_or_url.url if isinstance(ref_or_url, DatasetRef) else ref_or_url
    _throttle()
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    return _parse_dataset_file(resp.text)


def _parse_dataset_file(content: str) -> Dataset:
    """Parse the ``# header + numeric block`` exfortables_py format."""
    ds = Dataset()
    lines = content.splitlines()

    section: Optional[str] = None
    column_header: Optional[str] = None

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue

        if stripped.startswith("#"):
            payload = stripped.lstrip("#").strip()
            if not payload:
                continue

            # Section markers ("Header:", "Target:", ...)
            if payload.endswith(":") and ":" not in payload[:-1]:
                section = payload[:-1].strip().lower()
                continue

            # Header column line (last comment line before data)
            if any(tok in payload for tok in ("E_in", "E(MeV)", "Energy", "XS",
                                              "Angle", "Yield", "DATA")):
                if "(" in payload and ")" in payload:
                    column_header = payload

            # key: value within a section
            if ":" in payload:
                key, _, value = payload.partition(":")
                key = key.strip()
                value = value.strip()
                _assign_metadata(ds, section, key, value)
            continue

        # Numeric data row
        try:
            row = [float(x) for x in stripped.split()]
        except ValueError:
            continue
        if row:
            ds.rows.append(row)

    if column_header:
        # Each label is a contiguous run of non-whitespace characters
        # (e.g. "E_in(MeV)", "dE_in(MeV)", "XS(B)", "dXS(B)").
        ds.column_labels = column_header.split()

    return ds


def _assign_metadata(ds: Dataset, section: Optional[str], key: str, value: str) -> None:
    """Map header key/value pairs into Dataset fields."""
    k = key.lower()
    if section == "header":
        if k == "title":
            ds.title = value
    elif section == "target":
        if k == "z":
            ds.target_z = value
        elif k == "a":
            ds.target_a = value
        elif k == "nuclide":
            ds.target = value
    elif section == "reaction":
        if k == "process":
            ds.process = value
        elif "mf" in k and "mt" in k:
            # "MF-MT number: 3 - 18"
            ds.mf, _, ds.mt = value.partition("-")
            ds.mf = ds.mf.strip()
            ds.mt = ds.mt.strip()
        elif k.startswith("incident energy"):
            ds.energy_range = value
    elif section == "residual":
        if k == "nuclide":
            ds.residual = value
    elif section in ("exfor bib", "bib"):
        if k == "entry id":
            ds.entry_id = value.split()[0] if value else ""
        elif k == "reaction code":
            ds.reaction_code = value
        elif k.startswith("first author") or k == "author":
            ds.author = value
        elif k == "institute":
            ds.institute = value
        elif k == "reference":
            ds.reference = value
        elif k == "year":
            ds.year = value
        elif k == "facility":
            ds.facility = value
        elif k == "master file":
            ds.master_file_url = value
        elif k == "nds":
            ds.nds_url = value


# ── Full ENTRY JSON retrieval ─────────────────────────────────────────


def get_entry_json(entry_id: str) -> dict:
    """Fetch a complete EXFOR ENTRY as JSON from IAEA-NDS/exfor_json.

    Parameters
    ----------
    entry_id : str
        EXFOR entry accession number (e.g. ``"14015"`` or ``"21082"``).

    Returns
    -------
    dict
        Top-level keys::

            entry, last_updated, number_of_revisions, histories,
            bib_record, reactions, data_tables, experimental_conditions

        ``bib_record`` carries the full BIB block (TITLE, AUTHOR,
        INSTITUTE, FACILITY, REFERENCE, MONITOR, COMMENT, ...).
        ``data_tables`` carries every subentry's data array.

    Notes
    -----
    The repository organises files as ``json/{first_three_digits}/
    {entry_id}.json`` — this function builds that path automatically.
    """
    eid = entry_id.strip()
    if not eid:
        raise ValueError("entry_id must be non-empty")
    # Pad / take first 3 digits for the directory bucket
    bucket = eid[:3].zfill(3) if eid[:3].isdigit() else eid[:3]
    url = (
        f"{_GH_RAW}/{_EXFOR_JSON_OWNER}/{_EXFOR_JSON_REPO}/"
        f"main/json/{bucket}/{eid}.json"
    )
    _throttle()
    resp = requests.get(url, timeout=30)
    if resp.status_code == 404:
        raise LookupError(f"No EXFOR entry json/{bucket}/{eid}.json — entry "
                          f"may not exist or may not yet be in the JSON dump.")
    resp.raise_for_status()
    return resp.json()
