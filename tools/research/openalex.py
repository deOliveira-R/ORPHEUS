"""OpenAlex API client for searching scholarly works.

Provides typed access to the OpenAlex REST API for finding papers,
authors, and citation networks across 250M+ works. Includes open-access
PDF links, concept tagging, and full citation graphs.

No authentication required. Add email for polite pool (faster rate limits).

Usage::

    from tools.research.openalex import search, get_work, search_by_doi

    # Free-text search
    results = search("discrete ordinates neutron transport")

    # Filtered search
    results = search(
        "collision probability",
        publication_year=2020,
        open_access=True,
    )

    # Fetch by DOI
    work = get_work("https://doi.org/10.13182/NSE08-64")

    # Search with filters
    results = search(
        "Monte Carlo",
        publication_year="2015-2023",
        type="article",
    )
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Optional, Union

import requests

BASE_URL = "https://api.openalex.org"

# Optional: set to your email for polite pool (10x higher rate limit).
MAILTO: Optional[str] = None

_MIN_REQUEST_INTERVAL = 0.1  # 10 req/s for polite pool
_last_request_time = 0.0


def _throttle() -> None:
    global _last_request_time
    elapsed = time.monotonic() - _last_request_time
    if elapsed < _MIN_REQUEST_INTERVAL:
        time.sleep(_MIN_REQUEST_INTERVAL - elapsed)
    _last_request_time = time.monotonic()


def _params_base() -> dict:
    """Base parameters including polite pool email."""
    if MAILTO:
        return {"mailto": MAILTO}
    return {}


def _reconstruct_abstract(inverted_index: dict) -> str:
    """Reconstruct abstract text from OpenAlex inverted index format."""
    if not inverted_index:
        return ""
    word_positions: list[tuple[int, str]] = []
    for word, positions in inverted_index.items():
        for pos in positions:
            word_positions.append((pos, word))
    word_positions.sort()
    return " ".join(word for _, word in word_positions)


# ── Data classes ──────────────────────────────────────────────────────


@dataclass
class Author:
    """An author on an OpenAlex work."""

    name: str
    openalex_id: str = ""
    orcid: str = ""
    institution: str = ""


@dataclass
class Work:
    """A single OpenAlex scholarly work."""

    openalex_id: str
    title: str
    authors: list[Author] = field(default_factory=list)
    abstract: str = ""
    publication_year: int = 0
    publication_date: str = ""
    doi: str = ""
    cited_by_count: int = 0
    work_type: str = ""
    is_open_access: bool = False
    oa_url: str = ""
    journal: str = ""
    volume: str = ""
    issue: str = ""
    first_page: str = ""
    last_page: str = ""
    concepts: list[str] = field(default_factory=list)

    @classmethod
    def from_json(cls, data: dict) -> Work:
        """Parse a work from the OpenAlex API JSON response."""
        # Authors
        authors = []
        for authorship in data.get("authorships") or []:
            author_data = authorship.get("author") or {}
            institutions = authorship.get("institutions") or []
            inst_name = institutions[0].get("display_name", "") if institutions else ""
            authors.append(Author(
                name=author_data.get("display_name", ""),
                openalex_id=author_data.get("id", ""),
                orcid=author_data.get("orcid") or "",
                institution=inst_name,
            ))

        # Abstract
        abstract = _reconstruct_abstract(
            data.get("abstract_inverted_index") or {}
        )

        # Open access
        oa = data.get("open_access") or {}
        oa_url = oa.get("oa_url") or ""

        # Primary location (journal info)
        location = data.get("primary_location") or {}
        source = location.get("source") or {}
        journal = source.get("display_name", "")

        # Biblio
        biblio = data.get("biblio") or {}

        # Concepts
        concepts = []
        for concept in data.get("concepts") or []:
            if concept.get("score", 0) > 0.3:
                concepts.append(concept.get("display_name", ""))

        return cls(
            openalex_id=data.get("id", ""),
            title=data.get("title") or "",
            authors=authors,
            abstract=abstract,
            publication_year=data.get("publication_year") or 0,
            publication_date=data.get("publication_date") or "",
            doi=(data.get("doi") or "").replace("https://doi.org/", ""),
            cited_by_count=data.get("cited_by_count", 0),
            work_type=data.get("type", ""),
            is_open_access=oa.get("is_oa", False),
            oa_url=oa_url,
            journal=journal,
            volume=biblio.get("volume") or "",
            issue=biblio.get("issue") or "",
            first_page=biblio.get("first_page") or "",
            last_page=biblio.get("last_page") or "",
            concepts=concepts,
        )

    @property
    def first_author(self) -> str:
        return self.authors[0].name if self.authors else "Unknown"

    def summary(self) -> str:
        """One-line summary for display."""
        year = self.publication_year or "n.d."
        oa = " [OA]" if self.is_open_access else ""
        cite = f", cited {self.cited_by_count}x" if self.cited_by_count else ""
        return f"{self.first_author} ({year}): {self.title}{cite}{oa}"


@dataclass
class SearchResult:
    """Container for paginated OpenAlex search results."""

    works: list[Work]
    total_count: int
    page: int
    per_page: int

    @property
    def total_pages(self) -> int:
        if self.per_page <= 0:
            return 0
        return (self.total_count + self.per_page - 1) // self.per_page


# ── API functions ─────────────────────────────────────────────────────


def search(
    q: str,
    *,
    publication_year: Optional[Union[int, str]] = None,
    type: Optional[str] = None,
    open_access: Optional[bool] = None,
    per_page: int = 10,
    page: int = 1,
    sort: str = "relevance_score:desc",
) -> SearchResult:
    """Search OpenAlex for scholarly works.

    Parameters
    ----------
    q : str
        Search query (searches title, abstract, full text).
    publication_year : int or str, optional
        Single year (2020) or range ("2015-2023").
    type : str, optional
        Work type filter: "article", "review", "book-chapter",
        "proceedings-article", "dissertation", etc.
    open_access : bool, optional
        If True, only open-access works.
    per_page : int
        Results per page (default 10, max 200).
    page : int
        Page number (1-indexed).
    sort : str
        Sort: "relevance_score:desc", "cited_by_count:desc",
        "publication_date:desc", etc.

    Returns
    -------
    SearchResult
    """
    params = _params_base()
    params["search"] = q
    params["per_page"] = per_page
    params["page"] = page
    params["sort"] = sort

    # Build filter string
    filters = []
    if publication_year is not None:
        if isinstance(publication_year, int):
            filters.append(f"publication_year:{publication_year}")
        else:
            # Range: "2015-2023" → "publication_year:2015-2023"
            filters.append(f"publication_year:{publication_year}")
    if type is not None:
        filters.append(f"type:{type}")
    if open_access is True:
        filters.append("open_access.is_oa:true")
    if filters:
        params["filter"] = ",".join(filters)

    _throttle()
    resp = requests.get(f"{BASE_URL}/works", params=params, timeout=30)
    resp.raise_for_status()

    data = resp.json()
    meta = data.get("meta", {})
    total = meta.get("count", 0)

    works = [Work.from_json(item) for item in data.get("results", [])]

    return SearchResult(works=works, total_count=total, page=page, per_page=per_page)


def get_work(work_id: str) -> Work:
    """Fetch a single work by OpenAlex ID or DOI.

    Parameters
    ----------
    work_id : str
        OpenAlex ID (e.g. "W2741809807"), DOI URL
        (e.g. "https://doi.org/10.13182/NSE08-64"), or bare DOI.

    Returns
    -------
    Work
    """
    # Normalize bare DOI to URL form
    if work_id.startswith("10."):
        work_id = f"https://doi.org/{work_id}"

    _throttle()
    resp = requests.get(
        f"{BASE_URL}/works/{work_id}",
        params=_params_base(),
        timeout=30,
    )
    resp.raise_for_status()
    return Work.from_json(resp.json())


def get_citations(work_id: str, per_page: int = 25, page: int = 1) -> SearchResult:
    """Get works that cite a given work.

    Parameters
    ----------
    work_id : str
        OpenAlex ID or DOI of the cited work.
    per_page : int
        Results per page.
    page : int
        Page number.

    Returns
    -------
    SearchResult
    """
    if work_id.startswith("10."):
        work_id = f"https://doi.org/{work_id}"

    params = _params_base()
    params["filter"] = f"cites:{work_id}"
    params["per_page"] = per_page
    params["page"] = page
    params["sort"] = "cited_by_count:desc"

    _throttle()
    resp = requests.get(f"{BASE_URL}/works", params=params, timeout=30)
    resp.raise_for_status()

    data = resp.json()
    meta = data.get("meta", {})
    total = meta.get("count", 0)

    works = [Work.from_json(item) for item in data.get("results", [])]

    return SearchResult(works=works, total_count=total, page=page, per_page=per_page)


def get_references(work_id: str, per_page: int = 25, page: int = 1) -> SearchResult:
    """Get works referenced by a given work.

    Parameters
    ----------
    work_id : str
        OpenAlex ID or DOI of the citing work.
    per_page : int
        Results per page.
    page : int
        Page number.

    Returns
    -------
    SearchResult
    """
    if work_id.startswith("10."):
        work_id = f"https://doi.org/{work_id}"

    params = _params_base()
    params["filter"] = f"cited_by:{work_id}"
    params["per_page"] = per_page
    params["page"] = page

    _throttle()
    resp = requests.get(f"{BASE_URL}/works", params=params, timeout=30)
    resp.raise_for_status()

    data = resp.json()
    meta = data.get("meta", {})
    total = meta.get("count", 0)

    works = [Work.from_json(item) for item in data.get("results", [])]

    return SearchResult(works=works, total_count=total, page=page, per_page=per_page)
