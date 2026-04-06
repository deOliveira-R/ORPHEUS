"""CrossRef API client for DOI metadata lookup and scholarly work search.

Provides typed access to the CrossRef REST API for resolving DOIs,
searching publications, and querying journal-specific content.
Covers 150M+ scholarly works across all publishers.

No authentication required. Add email for polite pool (faster rate limits).

Usage::

    from tools.research.crossref import get_work, search, journal_search

    # Resolve a DOI
    work = get_work("10.13182/NSE08-64")

    # Search for papers
    results = search("collision probability neutron transport")

    # Search within a specific journal (by ISSN)
    results = journal_search("0029-5639", query="discrete ordinates")
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Optional

import requests

BASE_URL = "https://api.crossref.org"

# Optional: set to your email for polite pool.
MAILTO: Optional[str] = None

_MIN_REQUEST_INTERVAL = 0.1
_last_request_time = 0.0


def _throttle() -> None:
    global _last_request_time
    elapsed = time.monotonic() - _last_request_time
    if elapsed < _MIN_REQUEST_INTERVAL:
        time.sleep(_MIN_REQUEST_INTERVAL - elapsed)
    _last_request_time = time.monotonic()


def _headers() -> dict[str, str]:
    h = {"Accept": "application/json"}
    if MAILTO:
        h["User-Agent"] = f"ORPHEUS/1.0 (mailto:{MAILTO})"
    return h


def _params_base() -> dict:
    if MAILTO:
        return {"mailto": MAILTO}
    return {}


# ���─ Data classes ──────────────────────────────────────────────────────


@dataclass
class Author:
    """An author on a CrossRef work."""

    given: str = ""
    family: str = ""
    orcid: str = ""

    @property
    def name(self) -> str:
        return f"{self.given} {self.family}".strip()


@dataclass
class Work:
    """A single CrossRef work (journal article, book chapter, etc.)."""

    doi: str
    title: str
    authors: list[Author] = field(default_factory=list)
    abstract: str = ""
    container_title: str = ""  # journal name
    volume: str = ""
    issue: str = ""
    page: str = ""
    published_date: str = ""
    work_type: str = ""
    is_referenced_by_count: int = 0
    references_count: int = 0
    issn: list[str] = field(default_factory=list)
    publisher: str = ""
    url: str = ""
    license_url: str = ""

    @classmethod
    def from_json(cls, data: dict) -> Work:
        """Parse a work from the CrossRef API JSON response."""
        # Title
        titles = data.get("title") or []
        title = titles[0] if titles else ""

        # Authors
        authors = []
        for a in data.get("author") or []:
            authors.append(Author(
                given=a.get("given", ""),
                family=a.get("family", ""),
                orcid=a.get("ORCID", ""),
            ))

        # Abstract (may contain JATS XML tags)
        abstract = data.get("abstract", "")
        # Strip JATS XML tags if present
        if abstract and "<" in abstract:
            import re
            abstract = re.sub(r"<[^>]+>", "", abstract).strip()

        # Container (journal) title
        containers = data.get("container-title") or []
        container = containers[0] if containers else ""

        # Published date
        pub = data.get("published") or data.get("published-print") or {}
        date_parts = pub.get("date-parts", [[]])[0]
        if date_parts:
            published_date = "-".join(str(p).zfill(2) for p in date_parts)
        else:
            published_date = ""

        # License
        licenses = data.get("license") or []
        license_url = licenses[0].get("URL", "") if licenses else ""

        return cls(
            doi=data.get("DOI", ""),
            title=title,
            authors=authors,
            abstract=abstract,
            container_title=container,
            volume=data.get("volume", ""),
            issue=data.get("issue", ""),
            page=data.get("page", ""),
            published_date=published_date,
            work_type=data.get("type", ""),
            is_referenced_by_count=data.get("is-referenced-by-count", 0),
            references_count=data.get("references-count", 0),
            issn=data.get("ISSN") or [],
            publisher=data.get("publisher", ""),
            url=data.get("URL", ""),
            license_url=license_url,
        )

    @property
    def first_author(self) -> str:
        return self.authors[0].name if self.authors else "Unknown"

    def summary(self) -> str:
        """One-line summary for display."""
        year = self.published_date[:4] if self.published_date else "n.d."
        cite = f", cited {self.is_referenced_by_count}x" if self.is_referenced_by_count else ""
        return f"{self.first_author} ({year}): {self.title}{cite}"


@dataclass
class SearchResult:
    """Container for paginated CrossRef search results."""

    works: list[Work]
    total_count: int
    offset: int
    rows: int

    @property
    def total_pages(self) -> int:
        if self.rows <= 0:
            return 0
        return (self.total_count + self.rows - 1) // self.rows


# ── API functions ─────────────────────────────────────────────────────


def get_work(doi: str) -> Work:
    """Fetch a single work by DOI.

    Parameters
    ----------
    doi : str
        The DOI (e.g. "10.13182/NSE08-64").

    Returns
    -------
    Work
    """
    _throttle()
    resp = requests.get(
        f"{BASE_URL}/works/{doi}",
        headers=_headers(),
        params=_params_base(),
        timeout=30,
    )
    resp.raise_for_status()
    return Work.from_json(resp.json().get("message", {}))


def search(
    query: str,
    *,
    query_title: Optional[str] = None,
    query_author: Optional[str] = None,
    filter: Optional[str] = None,
    from_pub_date: Optional[str] = None,
    until_pub_date: Optional[str] = None,
    type: Optional[str] = None,
    has_abstract: Optional[bool] = None,
    rows: int = 10,
    offset: int = 0,
    sort: Optional[str] = None,
    order: str = "desc",
) -> SearchResult:
    """Search CrossRef for scholarly works.

    Parameters
    ----------
    query : str
        Free-text bibliographic search.
    query_title : str, optional
        Search within titles.
    query_author : str, optional
        Search by author name.
    filter : str, optional
        Raw CrossRef filter string (e.g. "from-pub-date:2020,type:journal-article").
    from_pub_date, until_pub_date : str, optional
        Date range in YYYY-MM-DD or YYYY format.
    type : str, optional
        Work type: "journal-article", "book-chapter", "proceedings-article", etc.
    has_abstract : bool, optional
        If True, only works with abstracts.
    rows : int
        Results per page (default 10, max 1000).
    offset : int
        Zero-based result offset.
    sort : str, optional
        Sort field: "relevance", "published", "is-referenced-by-count", etc.
    order : str
        Sort direction: "asc" or "desc".

    Returns
    -------
    SearchResult
    """
    params = _params_base()
    params["query"] = query
    params["rows"] = rows
    params["offset"] = offset
    if query_title:
        params["query.title"] = query_title
    if query_author:
        params["query.author"] = query_author
    if sort:
        params["sort"] = sort
        params["order"] = order

    # Build filter
    filter_parts = []
    if filter:
        filter_parts.append(filter)
    if from_pub_date:
        filter_parts.append(f"from-pub-date:{from_pub_date}")
    if until_pub_date:
        filter_parts.append(f"until-pub-date:{until_pub_date}")
    if type:
        filter_parts.append(f"type:{type}")
    if has_abstract is True:
        filter_parts.append("has-abstract:true")
    if filter_parts:
        params["filter"] = ",".join(filter_parts)

    _throttle()
    resp = requests.get(
        f"{BASE_URL}/works",
        headers=_headers(),
        params=params,
        timeout=30,
    )
    resp.raise_for_status()

    message = resp.json().get("message", {})
    total = message.get("total-results", 0)
    works = [Work.from_json(item) for item in message.get("items", [])]

    return SearchResult(works=works, total_count=total, offset=offset, rows=rows)


def journal_search(
    issn: str,
    query: str = "",
    *,
    rows: int = 10,
    offset: int = 0,
    sort: Optional[str] = None,
    order: str = "desc",
) -> SearchResult:
    """Search within a specific journal by ISSN.

    Parameters
    ----------
    issn : str
        Journal ISSN (e.g. "0029-5639" for Nuclear Science and Engineering).
    query : str
        Search query within the journal.
    rows : int
        Results per page.
    offset : int
        Zero-based result offset.
    sort : str, optional
        Sort field.
    order : str
        Sort direction.

    Returns
    -------
    SearchResult

    Common ISSNs for nuclear engineering journals::

        NSE  = "0029-5639"   # Nuclear Science and Engineering
        ANE  = "0306-4549"   # Annals of Nuclear Energy
        NED  = "0029-5493"   # Nuclear Engineering and Design
        NET  = "1738-5733"   # Nuclear Engineering and Technology
        JCP  = "0021-9991"   # Journal of Computational Physics
        PNE  = "0149-1970"   # Progress in Nuclear Energy
        JNM  = "0022-3115"   # Journal of Nuclear Materials
    """
    params = _params_base()
    if query:
        params["query"] = query
    params["rows"] = rows
    params["offset"] = offset
    if sort:
        params["sort"] = sort
        params["order"] = order

    _throttle()
    resp = requests.get(
        f"{BASE_URL}/journals/{issn}/works",
        headers=_headers(),
        params=params,
        timeout=30,
    )
    resp.raise_for_status()

    message = resp.json().get("message", {})
    total = message.get("total-results", 0)
    works = [Work.from_json(item) for item in message.get("items", [])]

    return SearchResult(works=works, total_count=total, offset=offset, rows=rows)
