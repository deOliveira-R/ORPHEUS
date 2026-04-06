"""Semantic Scholar API client for AI-enhanced paper search.

Provides typed access to the Semantic Scholar Academic Graph API
for finding papers, citation networks, and influential citations
across 200M+ papers.

No API key required for basic use (100 req/min). Free API key
available for higher throughput (https://www.semanticscholar.org/product/api).

Usage::

    from tools.research.semantic_scholar import search, get_paper, get_citations

    # Search for papers
    results = search("discrete ordinates neutron transport")

    # Fetch a paper by DOI, arXiv ID, or S2 ID
    paper = get_paper("DOI:10.13182/NSE08-64")
    paper = get_paper("ARXIV:2103.12345")
    paper = get_paper("649def34f8be52c8b66281af98ae884c09aef38b")

    # Get citing papers
    citations = get_citations("DOI:10.13182/NSE08-64")
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Optional

import requests

BASE_URL = "https://api.semanticscholar.org/graph/v1"

# Optional: set to your API key for higher rate limits.
API_KEY: Optional[str] = None

# Without key: 100 req/min. With key: 1000 req/min.
_MIN_REQUEST_INTERVAL = 0.6  # Conservative for no-key use
_last_request_time = 0.0

# Default fields to request for papers.
_PAPER_FIELDS = (
    "paperId,title,abstract,authors,year,citationCount,"
    "influentialCitationCount,isOpenAccess,openAccessPdf,"
    "externalIds,publicationTypes,journal,venue,referenceCount,"
    "fieldsOfStudy,s2FieldsOfStudy,publicationDate"
)

_AUTHOR_FIELDS = "authorId,name,affiliations,citationCount,hIndex"


def _throttle() -> None:
    global _last_request_time
    elapsed = time.monotonic() - _last_request_time
    if elapsed < _MIN_REQUEST_INTERVAL:
        time.sleep(_MIN_REQUEST_INTERVAL - elapsed)
    _last_request_time = time.monotonic()


def _headers() -> dict[str, str]:
    h: dict[str, str] = {}
    if API_KEY:
        h["x-api-key"] = API_KEY
    return h


# ── Data classes ──────────────────────────────────────────────────────


@dataclass
class Author:
    """An author on a Semantic Scholar paper."""

    name: str
    author_id: str = ""


@dataclass
class Paper:
    """A single Semantic Scholar paper."""

    paper_id: str
    title: str
    authors: list[Author] = field(default_factory=list)
    abstract: str = ""
    year: int = 0
    publication_date: str = ""
    citation_count: int = 0
    influential_citation_count: int = 0
    reference_count: int = 0
    is_open_access: bool = False
    oa_pdf_url: str = ""
    doi: str = ""
    arxiv_id: str = ""
    journal: str = ""
    venue: str = ""
    fields_of_study: list[str] = field(default_factory=list)
    publication_types: list[str] = field(default_factory=list)

    @classmethod
    def from_json(cls, data: dict) -> Paper:
        """Parse a paper from the S2 API JSON response."""
        authors = [
            Author(
                name=a.get("name", ""),
                author_id=a.get("authorId") or "",
            )
            for a in (data.get("authors") or [])
        ]

        ext_ids = data.get("externalIds") or {}
        oa_pdf = data.get("openAccessPdf") or {}

        journal_data = data.get("journal") or {}
        journal_name = journal_data.get("name", "") if isinstance(journal_data, dict) else ""

        fields = []
        for f in data.get("s2FieldsOfStudy") or data.get("fieldsOfStudy") or []:
            if isinstance(f, dict):
                fields.append(f.get("category", ""))
            elif isinstance(f, str):
                fields.append(f)

        return cls(
            paper_id=data.get("paperId", ""),
            title=data.get("title") or "",
            authors=authors,
            abstract=data.get("abstract") or "",
            year=data.get("year") or 0,
            publication_date=data.get("publicationDate") or "",
            citation_count=data.get("citationCount", 0),
            influential_citation_count=data.get("influentialCitationCount", 0),
            reference_count=data.get("referenceCount", 0),
            is_open_access=data.get("isOpenAccess", False),
            oa_pdf_url=oa_pdf.get("url", ""),
            doi=ext_ids.get("DOI", ""),
            arxiv_id=ext_ids.get("ArXiv", ""),
            journal=journal_name,
            venue=data.get("venue") or "",
            fields_of_study=fields,
            publication_types=data.get("publicationTypes") or [],
        )

    @property
    def first_author(self) -> str:
        return self.authors[0].name if self.authors else "Unknown"

    def summary(self) -> str:
        """One-line summary for display."""
        year = self.year or "n.d."
        cite = f", cited {self.citation_count}x" if self.citation_count else ""
        inf = f" ({self.influential_citation_count} influential)" if self.influential_citation_count else ""
        oa = " [OA]" if self.is_open_access else ""
        return f"{self.first_author} ({year}): {self.title}{cite}{inf}{oa}"

    def s2_url(self) -> str:
        """Semantic Scholar paper page URL."""
        return f"https://www.semanticscholar.org/paper/{self.paper_id}"


@dataclass
class SearchResult:
    """Container for paginated Semantic Scholar search results."""

    papers: list[Paper]
    total_count: int
    offset: int
    limit: int

    @property
    def total_pages(self) -> int:
        if self.limit <= 0:
            return 0
        return (self.total_count + self.limit - 1) // self.limit


# ── API functions ─────────────────────────────────────────────────────


def search(
    query: str,
    *,
    year: Optional[str] = None,
    fields_of_study: Optional[list[str]] = None,
    open_access_only: bool = False,
    limit: int = 10,
    offset: int = 0,
) -> SearchResult:
    """Search Semantic Scholar for papers.

    Parameters
    ----------
    query : str
        Search query (title + abstract).
    year : str, optional
        Year filter: "2020" (single), "2015-2023" (range),
        "2020-" (from), "-2015" (until).
    fields_of_study : list[str], optional
        Filter by field: "Physics", "Computer Science",
        "Engineering", "Mathematics", etc.
    open_access_only : bool
        If True, only open-access papers.
    limit : int
        Results per page (default 10, max 100).
    offset : int
        Zero-based result offset (max 1000).

    Returns
    -------
    SearchResult
    """
    params: dict = {
        "query": query,
        "limit": limit,
        "offset": offset,
        "fields": _PAPER_FIELDS,
    }
    if year:
        params["year"] = year
    if fields_of_study:
        params["fieldsOfStudy"] = ",".join(fields_of_study)
    if open_access_only:
        params["openAccessPdf"] = ""

    _throttle()
    resp = requests.get(
        f"{BASE_URL}/paper/search",
        headers=_headers(),
        params=params,
        timeout=30,
    )
    resp.raise_for_status()

    data = resp.json()
    total = data.get("total", 0)
    papers = [Paper.from_json(item) for item in data.get("data", [])]

    return SearchResult(papers=papers, total_count=total, offset=offset, limit=limit)


def get_paper(paper_id: str) -> Paper:
    """Fetch a single paper by Semantic Scholar ID, DOI, or arXiv ID.

    Parameters
    ----------
    paper_id : str
        Paper identifier. Formats:
        - S2 ID: "649def34f8be52c8b66281af98ae884c09aef38b"
        - DOI: "DOI:10.13182/NSE08-64"
        - arXiv: "ARXIV:2103.12345"
        - PubMed: "PMID:12345678"
        - Corpus ID: "CorpusID:12345"

    Returns
    -------
    Paper
    """
    _throttle()
    resp = requests.get(
        f"{BASE_URL}/paper/{paper_id}",
        headers=_headers(),
        params={"fields": _PAPER_FIELDS},
        timeout=30,
    )
    resp.raise_for_status()
    return Paper.from_json(resp.json())


def get_citations(
    paper_id: str,
    *,
    limit: int = 25,
    offset: int = 0,
) -> SearchResult:
    """Get papers that cite a given paper.

    Parameters
    ----------
    paper_id : str
        Paper identifier (S2 ID, DOI:..., ARXIV:...).
    limit : int
        Results per page (max 1000).
    offset : int
        Zero-based offset.

    Returns
    -------
    SearchResult
    """
    _throttle()
    resp = requests.get(
        f"{BASE_URL}/paper/{paper_id}/citations",
        headers=_headers(),
        params={
            "fields": _PAPER_FIELDS,
            "limit": limit,
            "offset": offset,
        },
        timeout=30,
    )
    resp.raise_for_status()

    data = resp.json()
    total = data.get("total", len(data.get("data", [])))
    papers = []
    for item in data.get("data", []):
        citing = item.get("citingPaper", {})
        if citing.get("paperId"):
            papers.append(Paper.from_json(citing))

    return SearchResult(papers=papers, total_count=total, offset=offset, limit=limit)


def get_references(
    paper_id: str,
    *,
    limit: int = 25,
    offset: int = 0,
) -> SearchResult:
    """Get papers referenced by a given paper.

    Parameters
    ----------
    paper_id : str
        Paper identifier (S2 ID, DOI:..., ARXIV:...).
    limit : int
        Results per page (max 1000).
    offset : int
        Zero-based offset.

    Returns
    -------
    SearchResult
    """
    _throttle()
    resp = requests.get(
        f"{BASE_URL}/paper/{paper_id}/references",
        headers=_headers(),
        params={
            "fields": _PAPER_FIELDS,
            "limit": limit,
            "offset": offset,
        },
        timeout=30,
    )
    resp.raise_for_status()

    data = resp.json()
    total = data.get("total", len(data.get("data", [])))
    papers = []
    for item in data.get("data", []):
        cited = item.get("citedPaper", {})
        if cited.get("paperId"):
            papers.append(Paper.from_json(cited))

    return SearchResult(papers=papers, total_count=total, offset=offset, limit=limit)
