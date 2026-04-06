"""Elsevier / Scopus API client for searching scientific publications.

Provides typed access to the Scopus Search API and Abstract Retrieval API
for finding journal articles, conference papers, and reviews.

Requires an API key from https://dev.elsevier.com/apikey/manage.
Store it in ``tools/research/elsevier.toml`` (gitignored). See ``elsevier.toml.example``.

Usage::

    from tools.research.elsevier import search, get_abstract

    # Search by topic
    results = search("TITLE-ABS-KEY(discrete ordinates) AND PUBYEAR > 2010")

    # Convenience search with keyword arguments
    results = search_simple(
        keywords="collision probability",
        author="Hébert",
        journal="Annals of Nuclear Energy",
    )

    # Fetch full abstract by DOI
    abstract = get_abstract(doi="10.1016/j.anucene.2023.01.001")

Scopus query syntax reference:

    Field codes:  TITLE, ABS, KEY, TITLE-ABS-KEY, AUTH, DOI, SRCTITLE,
                  PUBYEAR, DOCTYPE, OPENACCESS, AFFIL, ...
    Operators:    AND, OR, AND NOT
    Phrases:      {exact phrase} or "loose phrase"
    Wildcards:    * (multi-char), ? (single-char)
    Date filter:  PUBYEAR > 2010, PUBYEAR = 2023
    Doc types:    DOCTYPE(ar) = article, DOCTYPE(re) = review, DOCTYPE(cp) = conf paper
"""

from __future__ import annotations

import tomllib
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import requests

BASE_URL = "https://api.elsevier.com/content"

_CONFIG_PATH = Path(__file__).parent / "elsevier.toml"

# Rate limit: max 9 req/s → ~0.12s between requests. Use 0.15s for margin.
_MIN_REQUEST_INTERVAL = 0.15
_last_request_time = 0.0


def _throttle() -> None:
    """Enforce minimum interval between requests."""
    global _last_request_time
    elapsed = time.monotonic() - _last_request_time
    if elapsed < _MIN_REQUEST_INTERVAL:
        time.sleep(_MIN_REQUEST_INTERVAL - elapsed)
    _last_request_time = time.monotonic()


def _load_config() -> dict[str, str]:
    """Load API credentials from elsevier.toml."""
    if not _CONFIG_PATH.exists():
        raise FileNotFoundError(
            f"Elsevier config not found at {_CONFIG_PATH}.\n"
            f"Copy elsevier.toml.example to elsevier.toml and add your API key."
        )
    with open(_CONFIG_PATH, "rb") as f:
        config = tomllib.load(f)
    section = config.get("elsevier", {})
    api_key = section.get("api_key", "")
    if not api_key or api_key == "YOUR_API_KEY_HERE":
        raise ValueError("Set your Elsevier API key in tools/elsevier.toml")
    return section


def _headers() -> dict[str, str]:
    """Build request headers with API key."""
    config = _load_config()
    h = {
        "X-ELS-APIKey": config["api_key"],
        "Accept": "application/json",
    }
    if "inst_token" in config:
        h["X-ELS-Insttoken"] = config["inst_token"]
    return h


# ── Data classes ──────────────────────────────────────────────────────


@dataclass
class Author:
    """An author on a Scopus record."""

    name: str
    author_id: str = ""
    orcid: str = ""


@dataclass
class Record:
    """A Scopus search result entry."""

    scopus_id: str
    title: str
    authors: list[Author] = field(default_factory=list)
    description: str = ""
    publication_date: str = ""
    doi: str = ""
    journal: str = ""
    volume: str = ""
    issue: str = ""
    pages: str = ""
    doc_type: str = ""
    cited_by: int = 0
    keywords: str = ""
    open_access: bool = False
    scopus_url: str = ""

    @classmethod
    def from_entry(cls, entry: dict) -> Record:
        """Parse a record from a Scopus search JSON entry."""
        # Parse authors
        authors = []
        for a in entry.get("author") or []:
            authors.append(Author(
                name=a.get("authname", ""),
                author_id=a.get("authid", ""),
                orcid=a.get("orcid", ""),
            ))

        # Find Scopus URL from links
        scopus_url = ""
        for lnk in entry.get("link") or []:
            if lnk.get("@ref") == "scopus":
                scopus_url = lnk.get("@href", "")
                break

        return cls(
            scopus_id=entry.get("dc:identifier", "").replace("SCOPUS_ID:", ""),
            title=entry.get("dc:title", ""),
            authors=authors,
            description=entry.get("dc:description", ""),
            publication_date=entry.get("prism:coverDate", ""),
            doi=entry.get("prism:doi", ""),
            journal=entry.get("prism:publicationName", ""),
            volume=entry.get("prism:volume", ""),
            issue=entry.get("prism:issueIdentifier", ""),
            pages=entry.get("prism:pageRange", ""),
            doc_type=entry.get("subtypeDescription", ""),
            cited_by=int(entry.get("citedby-count", 0)),
            keywords=entry.get("authkeywords", ""),
            open_access=entry.get("openaccessFlag") == "true",
            scopus_url=scopus_url,
        )

    @property
    def first_author(self) -> str:
        return self.authors[0].name if self.authors else "Unknown"

    def summary(self) -> str:
        """One-line summary for display."""
        year = self.publication_date[:4] if self.publication_date else "n.d."
        cite = f", cited {self.cited_by}x" if self.cited_by else ""
        return f"{self.first_author} ({year}): {self.title}{cite}"


@dataclass
class Abstract:
    """Full abstract retrieved from the Abstract Retrieval API."""

    scopus_id: str
    title: str
    authors: list[Author] = field(default_factory=list)
    abstract: str = ""
    doi: str = ""
    journal: str = ""
    volume: str = ""
    issue: str = ""
    pages: str = ""
    publication_date: str = ""
    cited_by: int = 0
    keywords: list[str] = field(default_factory=list)
    subject_areas: list[str] = field(default_factory=list)

    @classmethod
    def from_response(cls, data: dict) -> Abstract:
        """Parse from abstracts-retrieval-response JSON."""
        core = data.get("coredata", {})

        # Authors
        authors = []
        authors_data = data.get("authors", {}).get("author") or []
        for a in authors_data:
            pref = a.get("preferred-name") or {}
            name = f"{pref.get('ce:given-name', '')} {pref.get('ce:surname', '')}".strip()
            if not name:
                name = a.get("ce:indexed-name", "")
            authors.append(Author(
                name=name,
                author_id=a.get("@auid", ""),
            ))

        # Keywords
        kw_data = data.get("authkeywords", {})
        keywords = []
        if kw_data:
            for kw in kw_data.get("author-keyword") or []:
                keywords.append(kw.get("$", ""))

        # Subject areas
        subj_data = data.get("subject-areas", {})
        subject_areas = []
        if subj_data:
            for s in subj_data.get("subject-area") or []:
                subject_areas.append(s.get("$", ""))

        return cls(
            scopus_id=core.get("dc:identifier", "").replace("SCOPUS_ID:", ""),
            title=core.get("dc:title", ""),
            authors=authors,
            abstract=core.get("dc:description", ""),
            doi=core.get("prism:doi", ""),
            journal=core.get("prism:publicationName", ""),
            volume=core.get("prism:volume", ""),
            issue=core.get("prism:issueIdentifier", ""),
            pages=core.get("prism:pageRange", ""),
            publication_date=core.get("prism:coverDate", ""),
            cited_by=int(core.get("citedby-count", 0)),
            keywords=keywords,
            subject_areas=subject_areas,
        )


@dataclass
class SearchResult:
    """Container for paginated Scopus search results."""

    records: list[Record]
    total_count: int
    start: int
    count: int

    @property
    def total_pages(self) -> int:
        if self.count <= 0:
            return 0
        return (self.total_count + self.count - 1) // self.count


@dataclass
class QuotaInfo:
    """API quota status from response headers."""

    limit: int
    remaining: int
    reset_epoch: int

    @classmethod
    def from_headers(cls, headers: dict) -> Optional[QuotaInfo]:
        try:
            return cls(
                limit=int(headers.get("X-RateLimit-Limit", 0)),
                remaining=int(headers.get("X-RateLimit-Remaining", 0)),
                reset_epoch=int(headers.get("X-RateLimit-Reset", 0)),
            )
        except (ValueError, TypeError):
            return None


# ── API functions ─────────────────────────────────────────────────────


def search(
    query: str,
    *,
    start: int = 0,
    count: int = 25,
    view: str = "STANDARD",
    sort: Optional[str] = None,
    date: Optional[str] = None,
    subject: Optional[str] = None,
) -> SearchResult:
    """Search Scopus with a raw query string.

    Parameters
    ----------
    query : str
        Scopus boolean query (e.g. "TITLE-ABS-KEY(neutron transport)").
    start : int
        Zero-based offset (default 0).
    count : int
        Results per page (default 25).
    view : str
        "STANDARD" or "COMPLETE" (includes abstract + full author list).
    sort : str, optional
        Sort field with +/- prefix (e.g. "-citedby-count", "+coverDate").
    date : str, optional
        Year range filter (e.g. "2015-2023").
    subject : str, optional
        Subject area code (e.g. "ENER", "PHYS", "ENGI").

    Returns
    -------
    SearchResult
    """
    params: dict = {
        "query": query,
        "start": start,
        "count": count,
        "view": view,
    }
    if sort is not None:
        params["sort"] = sort
    if date is not None:
        params["date"] = date
    if subject is not None:
        params["subj"] = subject

    _throttle()
    resp = requests.get(
        f"{BASE_URL}/search/scopus",
        headers=_headers(),
        params=params,
        timeout=30,
    )
    resp.raise_for_status()

    data = resp.json().get("search-results", {})
    total = int(data.get("opensearch:totalResults", 0))
    entries = data.get("entry", [])

    # Scopus returns an error entry when no results match
    records = []
    for entry in entries:
        if "error" not in entry:
            records.append(Record.from_entry(entry))

    return SearchResult(records=records, total_count=total, start=start, count=count)


def search_simple(
    keywords: Optional[str] = None,
    *,
    title: Optional[str] = None,
    author: Optional[str] = None,
    journal: Optional[str] = None,
    doi: Optional[str] = None,
    year_from: Optional[int] = None,
    year_to: Optional[int] = None,
    open_access: Optional[bool] = None,
    start: int = 0,
    count: int = 25,
    sort: Optional[str] = None,
) -> SearchResult:
    """Build and execute a Scopus query from keyword arguments.

    Parameters
    ----------
    keywords : str, optional
        Search in title + abstract + keywords.
    title : str, optional
        Search in title only.
    author : str, optional
        Search by author name.
    journal : str, optional
        Filter by journal/source title.
    doi : str, optional
        Search by DOI.
    year_from, year_to : int, optional
        Publication year range.
    open_access : bool, optional
        If True, only open access articles.
    start, count, sort
        Pagination and sorting (see :func:`search`).

    Returns
    -------
    SearchResult
    """
    parts: list[str] = []
    if keywords is not None:
        parts.append(f"TITLE-ABS-KEY({keywords})")
    if title is not None:
        parts.append(f"TITLE({title})")
    if author is not None:
        parts.append(f"AUTH({author})")
    if journal is not None:
        parts.append(f'SRCTITLE("{journal}")')
    if doi is not None:
        parts.append(f"DOI({doi})")
    if year_from is not None:
        parts.append(f"PUBYEAR > {year_from - 1}")
    if year_to is not None:
        parts.append(f"PUBYEAR < {year_to + 1}")
    if open_access is True:
        parts.append("OPENACCESS(1)")

    if not parts:
        raise ValueError("At least one search term is required")

    query = " AND ".join(parts)
    return search(query, start=start, count=count, sort=sort)


def get_abstract(
    *,
    doi: Optional[str] = None,
    scopus_id: Optional[str] = None,
    eid: Optional[str] = None,
) -> Abstract:
    """Retrieve a full abstract from Scopus.

    Provide exactly one identifier.

    Parameters
    ----------
    doi : str, optional
        DOI (e.g. "10.1016/j.anucene.2023.01.001").
    scopus_id : str, optional
        Scopus ID.
    eid : str, optional
        EID (e.g. "2-s2.0-85012345678").

    Returns
    -------
    Abstract
    """
    if doi is not None:
        url = f"{BASE_URL}/abstract/doi/{doi}"
    elif scopus_id is not None:
        url = f"{BASE_URL}/abstract/scopus_id/{scopus_id}"
    elif eid is not None:
        url = f"{BASE_URL}/abstract/eid/{eid}"
    else:
        raise ValueError("Provide one of: doi, scopus_id, eid")

    _throttle()
    resp = requests.get(
        url,
        headers=_headers(),
        params={"view": "META_ABS"},
        timeout=30,
    )
    resp.raise_for_status()

    data = resp.json().get("abstracts-retrieval-response", {})
    return Abstract.from_response(data)


def quota() -> Optional[QuotaInfo]:
    """Check remaining API quota with a minimal request.

    Returns
    -------
    QuotaInfo or None
        Quota limits and remaining count, or None if unavailable.
    """
    _throttle()
    resp = requests.get(
        f"{BASE_URL}/search/scopus",
        headers=_headers(),
        params={"query": "TITLE(test)", "count": 1, "view": "STANDARD"},
        timeout=30,
    )
    resp.raise_for_status()
    return QuotaInfo.from_headers(resp.headers)
