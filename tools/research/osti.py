"""OSTI.gov API client for searching DOE scientific publications.

Provides typed access to the OSTI.gov REST API (v1) for finding
technical reports, journal articles, and other DOE-sponsored literature.

No authentication required. All endpoints are public GET requests.

Usage::

    from tools.research.osti import search, get_record, get_bibtex

    # Free-text search
    results = search("discrete ordinates cylindrical geometry")

    # Filtered search
    results = search(
        title="collision probability",
        author="Hébert",
        rows=5,
    )

    # Fetch a specific record
    record = get_record("1893820")

    # Get BibTeX for citation
    bibtex = get_bibtex("1893820")
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Optional
from urllib.parse import urljoin

import requests

BASE_URL = "https://www.osti.gov/api/v1/records"

# Polite rate limit: minimum seconds between requests.
_MIN_REQUEST_INTERVAL = 0.5
_last_request_time = 0.0


def _throttle() -> None:
    """Enforce minimum interval between requests."""
    global _last_request_time
    elapsed = time.monotonic() - _last_request_time
    if elapsed < _MIN_REQUEST_INTERVAL:
        time.sleep(_MIN_REQUEST_INTERVAL - elapsed)
    _last_request_time = time.monotonic()


# ── Data classes ──────────────────────────────────────────────────────


@dataclass
class Link:
    """A link associated with an OSTI record (citation page, fulltext PDF)."""

    rel: str
    href: str


@dataclass
class Record:
    """A single OSTI bibliographic record."""

    osti_id: str
    title: str
    authors: list[str] = field(default_factory=list)
    description: str = ""
    publication_date: str = ""
    doi: str = ""
    journal_name: str = ""
    journal_volume: str = ""
    journal_issue: str = ""
    product_type: str = ""
    report_number: str = ""
    research_orgs: list[str] = field(default_factory=list)
    sponsor_orgs: list[str] = field(default_factory=list)
    subjects: list[str] = field(default_factory=list)
    language: str = ""
    country: str = ""
    links: list[Link] = field(default_factory=list)

    @classmethod
    def from_json(cls, data: dict) -> Record:
        """Parse a record from the OSTI JSON response."""
        links = [
            Link(rel=lnk.get("rel", ""), href=lnk.get("href", ""))
            for lnk in (data.get("links") or [])
        ]
        return cls(
            osti_id=str(data.get("osti_id", "")),
            title=data.get("title", ""),
            authors=data.get("authors") or [],
            description=data.get("description", ""),
            publication_date=data.get("publication_date", ""),
            doi=data.get("doi", ""),
            journal_name=data.get("journal_name", ""),
            journal_volume=data.get("journal_volume", ""),
            journal_issue=data.get("journal_issue", ""),
            product_type=data.get("product_type", ""),
            report_number=data.get("report_number", ""),
            research_orgs=data.get("research_orgs") or [],
            sponsor_orgs=data.get("sponsor_orgs") or [],
            subjects=data.get("subjects") or [],
            language=data.get("language", ""),
            country=data.get("country_publication", ""),
            links=links,
        )

    @property
    def citation_url(self) -> Optional[str]:
        """URL of the OSTI citation page."""
        for lnk in self.links:
            if lnk.rel == "citation":
                return lnk.href
        return None

    @property
    def fulltext_url(self) -> Optional[str]:
        """URL of the fulltext PDF, if available."""
        for lnk in self.links:
            if lnk.rel == "fulltext":
                return lnk.href
        return None

    def summary(self) -> str:
        """One-line summary for display."""
        year = self.publication_date[:4] if self.publication_date else "n.d."
        first_author = self.authors[0].split(",")[0] if self.authors else "Unknown"
        return f"[{self.osti_id}] {first_author} ({year}): {self.title}"


@dataclass
class SearchResult:
    """Container for paginated search results."""

    records: list[Record]
    total_count: int
    page: int
    rows: int

    @property
    def total_pages(self) -> int:
        if self.rows <= 0:
            return 0
        return (self.total_count + self.rows - 1) // self.rows


# ── API functions ─────────────────────────────────────────────────────


def search(
    q: Optional[str] = None,
    *,
    title: Optional[str] = None,
    author: Optional[str] = None,
    subject: Optional[str] = None,
    fulltext: Optional[str] = None,
    biblio: Optional[str] = None,
    doi: Optional[str] = None,
    publication_date_start: Optional[str] = None,
    publication_date_end: Optional[str] = None,
    has_fulltext: Optional[bool] = None,
    rows: int = 20,
    page: int = 1,
    sort: str = "publication_date",
    order: str = "desc",
) -> SearchResult:
    """Search OSTI.gov for publications.

    Parameters
    ----------
    q : str, optional
        Free-text search across all fields.
    title, author, subject, fulltext, biblio : str, optional
        Field-specific searches.
    doi : str, optional
        Filter by DOI.
    publication_date_start, publication_date_end : str, optional
        Date range filter in MM/DD/YYYY format.
    has_fulltext : bool, optional
        If True, only return records with fulltext available.
    rows : int
        Results per page (default 20, max ~100).
    page : int
        Page number (1-indexed).
    sort : str
        Sort field (default "publication_date").
    order : str
        Sort direction: "asc" or "desc" (default "desc").

    Returns
    -------
    SearchResult
        Paginated results with total count.
    """
    params: dict = {}
    if q is not None:
        params["q"] = q
    if title is not None:
        params["title"] = title
    if author is not None:
        params["author"] = author
    if subject is not None:
        params["subject"] = subject
    if fulltext is not None:
        params["fulltext"] = fulltext
    if biblio is not None:
        params["biblio"] = biblio
    if doi is not None:
        params["doi"] = doi
    if publication_date_start is not None:
        params["publication_date_start"] = publication_date_start
    if publication_date_end is not None:
        params["publication_date_end"] = publication_date_end
    if has_fulltext is not None:
        params["has_fulltext"] = str(has_fulltext).lower()
    params["rows"] = rows
    params["page"] = page
    params["sort"] = sort
    params["order"] = order

    _throttle()
    resp = requests.get(
        BASE_URL, params=params, headers={"Accept": "application/json"}, timeout=30
    )
    resp.raise_for_status()

    total = int(resp.headers.get("X-Total-Count", 0))
    records = [Record.from_json(item) for item in resp.json()]

    return SearchResult(records=records, total_count=total, page=page, rows=rows)


def get_record(osti_id: str) -> Record:
    """Fetch a single OSTI record by ID.

    Parameters
    ----------
    osti_id : str
        The OSTI unique identifier.

    Returns
    -------
    Record
    """
    _throttle()
    resp = requests.get(
        f"{BASE_URL}/{osti_id}",
        headers={"Accept": "application/json"},
        timeout=30,
    )
    resp.raise_for_status()
    data = resp.json()
    if not data:
        raise ValueError(f"No record found for OSTI ID {osti_id}")
    return Record.from_json(data[0])


def get_bibtex(osti_id: str) -> str:
    """Fetch a BibTeX citation for an OSTI record.

    Parameters
    ----------
    osti_id : str
        The OSTI unique identifier.

    Returns
    -------
    str
        BibTeX entry as a string.
    """
    _throttle()
    resp = requests.get(
        f"{BASE_URL}/{osti_id}",
        headers={"Accept": "application/x-bibtex"},
        timeout=30,
    )
    resp.raise_for_status()
    return resp.text
