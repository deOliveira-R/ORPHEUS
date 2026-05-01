"""Zenodo REST API client.

Provides typed access to the CERN-hosted Zenodo records API. Zenodo
is the general-purpose research repository operated by CERN under
the EU OpenAIRE programme. Every deposit gets a citable DOI, so
Zenodo is the canonical home for:

- **Open-source nuclear code releases** with version DOIs:
  OpenMC, Serpent auxiliaries, DRAGON5/PARTISN inputs, NJOY/AMPX
  outputs, neutronics-related Python packages.
- **Cross-section / nuclear data library snapshots** that have a
  DOI but live nowhere else (e.g. JEFF processed libraries
  released by NEA Data Bank, ad-hoc covariance matrices).
- **Conference proceedings** (PHYSOR, M&C, ICAPP, NURETH, ICONE)
  that authors deposit alongside the published paper.
- **NEA / OECD benchmark input decks** (ICSBEP, IRPhE, BEAVRS
  variants) when a contributor uploads them.
- **Datasets** linked to a paper — measured neutron spectra,
  uncertainty matrices, fission-yield evaluations.

Zenodo runs InvenioRDM (same software family as INIS), so the API
shape is similar but the search syntax is Elasticsearch query string,
not Solr.

No authentication required for read operations. Add an API token
(``ZENODO_TOKEN`` env var) to lift the rate limit from 60/min to
100/min.

Usage::

    from tools.research.zenodo import search, get_record, get_bibtex

    # Free-text search
    results = search("OpenMC neutron transport")

    # Filter by resource type and Zenodo community
    results = search("MOC pin power", resource_type="dataset")
    results = search("nuclear data", communities="openmc")

    # Sort by recency
    results = search("cross section uncertainty", sort="mostrecent")

    # Resolve by Zenodo record id or DOI
    rec = get_record(12345)
    rec = get_record(doi="10.5281/zenodo.7233164")

    # BibTeX export
    print(get_bibtex(12345))
"""

from __future__ import annotations

import os
import time
from dataclasses import dataclass, field
from typing import Optional, Union

import requests

BASE_URL = "https://zenodo.org/api/records"
TOKEN: Optional[str] = os.environ.get("ZENODO_TOKEN")

# Polite throttle: stay well under the 60 req/min anonymous limit.
_MIN_REQUEST_INTERVAL = 1.1
_last_request_time = 0.0


def _throttle() -> None:
    global _last_request_time
    elapsed = time.monotonic() - _last_request_time
    if elapsed < _MIN_REQUEST_INTERVAL:
        time.sleep(_MIN_REQUEST_INTERVAL - elapsed)
    _last_request_time = time.monotonic()


def _auth_params() -> dict:
    return {"access_token": TOKEN} if TOKEN else {}


# ── Data classes ──────────────────────────────────────────────────────


@dataclass
class Creator:
    """A Zenodo record creator (author / depositor)."""

    name: str
    affiliation: str = ""
    orcid: str = ""


@dataclass
class FileEntry:
    """A file attached to a Zenodo record."""

    key: str  # filename
    size: int = 0
    checksum: str = ""
    download_url: str = ""


@dataclass
class Record:
    """A single Zenodo record (publication, dataset, software, ...)."""

    record_id: int
    doi: str = ""
    concept_doi: str = ""  # version-independent DOI
    title: str = ""
    creators: list[Creator] = field(default_factory=list)
    description: str = ""  # often HTML
    publication_date: str = ""
    resource_type: str = ""  # e.g. "publication", "software", "dataset"
    resource_subtype: str = ""  # e.g. "article", "report", "thesis"
    keywords: list[str] = field(default_factory=list)
    communities: list[str] = field(default_factory=list)
    license: str = ""
    version: str = ""
    is_open_access: bool = False
    html_url: str = ""
    files: list[FileEntry] = field(default_factory=list)
    related_identifiers: list[dict] = field(default_factory=list)
    journal_title: str = ""  # for resource_type=publication/article
    journal_volume: str = ""
    journal_issue: str = ""
    journal_pages: str = ""
    conference_title: str = ""
    notes: str = ""

    @classmethod
    def from_json(cls, data: dict) -> Record:
        """Parse a record from the Zenodo records API JSON response."""
        m = data.get("metadata", {}) or {}

        creators = []
        for c in m.get("creators", []) or []:
            creators.append(Creator(
                name=c.get("name", ""),
                affiliation=c.get("affiliation", "") or "",
                orcid=c.get("orcid", "") or "",
            ))

        files = []
        for f in data.get("files", []) or []:
            links = f.get("links", {}) or {}
            files.append(FileEntry(
                key=f.get("key", ""),
                size=f.get("size", 0) or 0,
                checksum=f.get("checksum", "") or "",
                download_url=links.get("self", "") or "",
            ))

        rt = m.get("resource_type", {}) or {}
        comm = [c.get("id", "") for c in (m.get("communities") or []) if c.get("id")]
        lic = (m.get("license") or {})
        if isinstance(lic, dict):
            lic_str = lic.get("id", "") or lic.get("title", "") or ""
        else:
            lic_str = str(lic)

        access = data.get("access", {}) or {}
        is_oa = access.get("status") == "open" or m.get("access_right") == "open"

        journal = m.get("journal", {}) or {}
        conference = m.get("conference", {}) or {}

        return cls(
            record_id=int(data.get("id", 0) or 0),
            doi=data.get("doi", "") or m.get("doi", "") or "",
            concept_doi=data.get("conceptdoi", "") or "",
            title=m.get("title", "") or "",
            creators=creators,
            description=m.get("description", "") or "",
            publication_date=m.get("publication_date", "") or "",
            resource_type=rt.get("type", "") or "",
            resource_subtype=rt.get("subtype", "") or "",
            keywords=list(m.get("keywords") or []),
            communities=comm,
            license=lic_str,
            version=m.get("version", "") or "",
            is_open_access=bool(is_oa),
            html_url=(data.get("links") or {}).get("self_html", "") or "",
            files=files,
            related_identifiers=list(m.get("related_identifiers") or []),
            journal_title=journal.get("title", "") or "",
            journal_volume=journal.get("volume", "") or "",
            journal_issue=journal.get("issue", "") or "",
            journal_pages=journal.get("pages", "") or "",
            conference_title=conference.get("title", "") or "",
            notes=m.get("notes", "") or "",
        )

    @property
    def first_creator(self) -> str:
        return self.creators[0].name if self.creators else "Unknown"

    @property
    def publication_year(self) -> str:
        return self.publication_date[:4] if self.publication_date else ""

    def summary(self) -> str:
        """One-line summary for display."""
        year = self.publication_year or "n.d."
        rt = f" [{self.resource_type}" + (f"/{self.resource_subtype}" if self.resource_subtype else "") + "]" if self.resource_type else ""
        oa = " [OA]" if self.is_open_access else ""
        return f"[Zenodo:{self.record_id}] {self.first_creator} ({year}){rt}: {self.title}{oa}"


@dataclass
class SearchResult:
    """Container for paginated Zenodo search results."""

    records: list[Record]
    total_count: int
    page: int
    size: int

    @property
    def total_pages(self) -> int:
        if self.size <= 0:
            return 0
        return (self.total_count + self.size - 1) // self.size


# ── API functions ─────────────────────────────────────────────────────


def search(
    q: str = "",
    *,
    resource_type: Optional[str] = None,
    communities: Optional[Union[str, list[str]]] = None,
    sort: str = "bestmatch",
    size: int = 10,
    page: int = 1,
    all_versions: bool = False,
) -> SearchResult:
    """Search Zenodo for records.

    Parameters
    ----------
    q : str
        Elasticsearch query string. Plain text is a free-text search;
        field-prefixed terms target Zenodo metadata fields. Useful
        prefixes::

            title:"OpenMC"                  exact title match
            creators.name:Romano            author surname
            keywords:"cross section"        keyword
            doi:10.5281/zenodo.7233164      DOI
            description:"Monte Carlo"       abstract/description text
            publication_date:[2020 TO 2024] date range
            resource_type.type:software     filter by type

    resource_type : str, optional
        Convenience filter on ``metadata.resource_type.type``. Common
        values: ``"publication"``, ``"poster"``, ``"presentation"``,
        ``"dataset"``, ``"software"``, ``"image"``, ``"video"``,
        ``"lesson"``, ``"physicalobject"``, ``"workflow"``, ``"other"``.
    communities : str or list[str], optional
        Restrict to one or more Zenodo community identifiers.
        Examples: ``"openmc"``, ``"oecd-nea-data-bank"``, ``"physor"``.
    sort : str
        ``"bestmatch"`` (relevance) or ``"mostrecent"``. Prefix with
        ``-`` for descending: ``"-mostrecent"``.
    size : int
        Results per page (max 25 anonymous, 100 with token).
    page : int
        Page number (1-indexed).
    all_versions : bool
        If True, return every version of each record (otherwise only
        the latest version of each concept DOI).

    Returns
    -------
    SearchResult
    """
    params = _auth_params()
    if q:
        params["q"] = q
    params["sort"] = sort
    params["size"] = size
    params["page"] = page
    if all_versions:
        params["all_versions"] = "true"
    if resource_type:
        params["type"] = resource_type
    if communities:
        comms = [communities] if isinstance(communities, str) else list(communities)
        params["communities"] = ",".join(comms)

    _throttle()
    resp = requests.get(BASE_URL, params=params, timeout=30)
    resp.raise_for_status()

    data = resp.json()
    hits = data.get("hits", {}) or {}
    total = hits.get("total", 0) or 0
    if isinstance(total, dict):
        total = total.get("value", 0)

    records = [Record.from_json(item) for item in hits.get("hits", []) or []]
    return SearchResult(records=records, total_count=total, page=page, size=size)


def get_record(
    record_id: Optional[Union[int, str]] = None,
    *,
    doi: Optional[str] = None,
) -> Record:
    """Fetch a single Zenodo record by id or DOI.

    Parameters
    ----------
    record_id : int or str, optional
        Numeric Zenodo record id (e.g. 7233164). Strings are accepted.
    doi : str, optional
        Zenodo DOI (e.g. "10.5281/zenodo.7233164" or full URL form).
        For non-Zenodo DOIs (10.x where x ≠ 5281/zenodo) raises
        ValueError — use CrossRef instead.

    Returns
    -------
    Record

    Raises
    ------
    ValueError
        If neither (or both) of record_id / doi is provided, or if
        a non-Zenodo DOI is given.
    LookupError
        If no record matches.
    """
    if (record_id is None) == (doi is None):
        raise ValueError("Provide exactly one of record_id or doi.")

    if record_id is not None:
        rid = str(record_id).strip()
    else:
        d = doi.replace("https://doi.org/", "").replace("http://doi.org/", "").strip()
        prefix = "10.5281/zenodo."
        if not d.lower().startswith(prefix):
            raise ValueError(
                f"DOI {doi!r} is not a Zenodo DOI (expected prefix {prefix!r})."
            )
        rid = d[len(prefix):]

    _throttle()
    resp = requests.get(
        f"{BASE_URL}/{rid}", params=_auth_params(), timeout=30
    )
    if resp.status_code == 404:
        raise LookupError(f"No Zenodo record with id={rid}")
    resp.raise_for_status()
    return Record.from_json(resp.json())


def get_bibtex(record_id: Union[int, str]) -> str:
    """Fetch a BibTeX entry for a Zenodo record.

    Uses content negotiation on the records endpoint (``Accept:
    application/x-bibtex``).

    Parameters
    ----------
    record_id : int or str
        Numeric Zenodo record id.

    Returns
    -------
    str
        BibTeX entry text.
    """
    _throttle()
    resp = requests.get(
        f"{BASE_URL}/{record_id}",
        params=_auth_params(),
        headers={"Accept": "application/x-bibtex"},
        timeout=30,
    )
    resp.raise_for_status()
    return resp.text


def get_versions(record_id: Union[int, str], size: int = 25) -> SearchResult:
    """List every published version of a record's concept DOI.

    Useful for software releases — call ``get_record(...)`` first to
    retrieve any one version, then use this to enumerate the version
    history.

    Parameters
    ----------
    record_id : int or str
        Any record id from the version family.
    size : int
        Max versions to return.

    Returns
    -------
    SearchResult
        Records in publication-date order (newest first).
    """
    rec = get_record(record_id)
    if not rec.concept_doi:
        # Single-version record — return as a 1-element result.
        return SearchResult(records=[rec], total_count=1, page=1, size=size)
    # The conceptdoi field is a Zenodo DOI; extract its id.
    concept_id = rec.concept_doi.split("zenodo.")[-1]
    return search(
        q=f"conceptrecid:{concept_id}",
        sort="-mostrecent",
        size=size,
        all_versions=True,
    )
