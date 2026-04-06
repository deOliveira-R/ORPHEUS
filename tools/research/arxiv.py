"""arXiv API client for searching preprints and published papers.

Provides typed access to the arXiv query API (Atom feed) for finding
papers by title, author, abstract, or category.

No authentication required. Rate limit: 3 seconds between requests.

Usage::

    from tools.research.arxiv import search, get_paper

    # Free-text search
    results = search("discrete ordinates cylindrical geometry")

    # Field-specific search
    results = search(title="collision probability", author="Hébert", max_results=5)

    # Category-filtered search
    results = search("neutron transport", category="nucl-th")

    # Fetch specific papers by arXiv ID
    papers = get_papers("2103.12345", "hep-ex/0307015")

arXiv query syntax (for the raw `query` parameter):

    Field prefixes: ti:, au:, abs:, co:, jr:, cat:, rn:, all:
    Boolean ops:    AND, OR, ANDNOT (uppercase)
    Phrases:        ti:"collision probability"
    Example:        au:morel AND ti:discrete ordinates
"""

from __future__ import annotations

import re
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import Optional

import requests

BASE_URL = "https://export.arxiv.org/api/query"

_NS = {
    "atom": "http://www.w3.org/2005/Atom",
    "opensearch": "http://a9.com/-/spec/opensearch/1.1/",
    "arxiv": "http://arxiv.org/schemas/atom",
}

# arXiv requires >= 3s between requests.
_MIN_REQUEST_INTERVAL = 3.0
_last_request_time = 0.0


def _throttle() -> None:
    """Enforce minimum interval between requests."""
    global _last_request_time
    elapsed = time.monotonic() - _last_request_time
    if elapsed < _MIN_REQUEST_INTERVAL:
        time.sleep(_MIN_REQUEST_INTERVAL - elapsed)
    _last_request_time = time.monotonic()


def _clean_text(text: Optional[str]) -> str:
    """Collapse whitespace in XML text content."""
    if text is None:
        return ""
    return re.sub(r"\s+", " ", text).strip()


# ── Data classes ──────────────────────────────────────────────────────


@dataclass
class Paper:
    """A single arXiv paper."""

    arxiv_id: str
    title: str
    authors: list[str] = field(default_factory=list)
    summary: str = ""
    published: str = ""
    updated: str = ""
    primary_category: str = ""
    categories: list[str] = field(default_factory=list)
    doi: str = ""
    journal_ref: str = ""
    comment: str = ""
    pdf_url: str = ""
    abstract_url: str = ""

    @classmethod
    def from_entry(cls, entry: ET.Element) -> Paper:
        """Parse a paper from an Atom <entry> element."""
        raw_id = entry.findtext("atom:id", default="", namespaces=_NS)
        # Extract arXiv ID from URL: http://arxiv.org/abs/2103.12345v1
        arxiv_id = raw_id.rsplit("/abs/", 1)[-1] if "/abs/" in raw_id else raw_id

        authors = [
            _clean_text(a.findtext("atom:name", namespaces=_NS))
            for a in entry.findall("atom:author", _NS)
        ]

        categories = [
            c.get("term", "")
            for c in entry.findall("atom:category", _NS)
        ]

        primary = entry.find("arxiv:primary_category", _NS)
        primary_cat = primary.get("term", "") if primary is not None else ""

        # Extract links by rel/title
        pdf_url = ""
        abstract_url = ""
        for link in entry.findall("atom:link", _NS):
            rel = link.get("rel", "")
            title = link.get("title", "")
            href = link.get("href", "")
            if title == "pdf":
                pdf_url = href
            elif rel == "alternate":
                abstract_url = href

        return cls(
            arxiv_id=arxiv_id,
            title=_clean_text(entry.findtext("atom:title", namespaces=_NS)),
            authors=authors,
            summary=_clean_text(entry.findtext("atom:summary", namespaces=_NS)),
            published=entry.findtext("atom:published", default="", namespaces=_NS),
            updated=entry.findtext("atom:updated", default="", namespaces=_NS),
            primary_category=primary_cat,
            categories=categories,
            doi=entry.findtext("arxiv:doi", default="", namespaces=_NS) or "",
            journal_ref=_clean_text(
                entry.findtext("arxiv:journal_ref", namespaces=_NS)
            ),
            comment=_clean_text(entry.findtext("arxiv:comment", namespaces=_NS)),
            pdf_url=pdf_url,
            abstract_url=abstract_url,
        )

    def display_id(self) -> str:
        """arXiv ID without version suffix."""
        return re.sub(r"v\d+$", "", self.arxiv_id)

    def summary_line(self) -> str:
        """One-line summary for display."""
        year = self.published[:4] if self.published else "n.d."
        first_author = self.authors[0].split(",")[0] if self.authors else "Unknown"
        return f"[{self.display_id()}] {first_author} ({year}): {self.title}"

    def bibtex_key(self) -> str:
        """Generate a BibTeX key from first author and year."""
        year = self.published[:4] if self.published else "nd"
        if self.authors:
            last = self.authors[0].split()[-1].lower()
            last = re.sub(r"[^a-z]", "", last)
        else:
            last = "unknown"
        return f"{last}{year}"


@dataclass
class SearchResult:
    """Container for paginated arXiv search results."""

    papers: list[Paper]
    total_count: int
    start: int
    max_results: int

    @property
    def total_pages(self) -> int:
        if self.max_results <= 0:
            return 0
        return (self.total_count + self.max_results - 1) // self.max_results


# ── Query builder ─────────────────────────────────────────────────────


def _build_query(
    q: Optional[str] = None,
    title: Optional[str] = None,
    author: Optional[str] = None,
    abstract: Optional[str] = None,
    category: Optional[str] = None,
) -> str:
    """Build an arXiv search_query string from individual fields."""
    parts: list[str] = []
    if q is not None:
        parts.append(f"all:{q}")
    if title is not None:
        parts.append(f"ti:{title}")
    if author is not None:
        parts.append(f"au:{author}")
    if abstract is not None:
        parts.append(f"abs:{abstract}")
    if category is not None:
        parts.append(f"cat:{category}")
    return " AND ".join(parts)


# ── Atom XML parser ──────────────────────────────────────────────────


def _parse_feed(xml_text: str) -> tuple[list[Paper], int]:
    """Parse an Atom feed into Papers and total result count."""
    root = ET.fromstring(xml_text)

    total = int(root.findtext("opensearch:totalResults", default="0", namespaces=_NS))

    papers = []
    for entry in root.findall("atom:entry", _NS):
        # arXiv returns a stub entry on errors — skip entries with no ID or "Error" title.
        entry_id = entry.findtext("atom:id", default="", namespaces=_NS)
        entry_title = entry.findtext("atom:title", default="", namespaces=_NS)
        if not entry_id or entry_title.strip() == "Error":
            continue
        papers.append(Paper.from_entry(entry))

    return papers, total


# ── API functions ─────────────────────────────────────────────────────


def search(
    q: Optional[str] = None,
    *,
    query: Optional[str] = None,
    title: Optional[str] = None,
    author: Optional[str] = None,
    abstract: Optional[str] = None,
    category: Optional[str] = None,
    start: int = 0,
    max_results: int = 10,
    sort_by: str = "relevance",
    sort_order: str = "descending",
) -> SearchResult:
    """Search arXiv for papers.

    Parameters
    ----------
    q : str, optional
        Free-text search across all fields.
    query : str, optional
        Raw arXiv query string (e.g. "au:morel AND ti:discrete ordinates").
        If provided, overrides q/title/author/abstract/category.
    title, author, abstract : str, optional
        Field-specific searches (combined with AND).
    category : str, optional
        arXiv category filter (e.g. "nucl-th", "physics.comp-ph").
    start : int
        Zero-based result offset (default 0).
    max_results : int
        Results per page (default 10, max 2000).
    sort_by : str
        "relevance", "lastUpdatedDate", or "submittedDate".
    sort_order : str
        "ascending" or "descending".

    Returns
    -------
    SearchResult
        Papers with total count and pagination info.
    """
    if query is not None:
        search_query = query
    else:
        search_query = _build_query(q, title, author, abstract, category)

    if not search_query:
        raise ValueError("At least one search term is required")

    params = {
        "search_query": search_query,
        "start": start,
        "max_results": max_results,
        "sortBy": sort_by,
        "sortOrder": sort_order,
    }

    _throttle()
    resp = requests.get(BASE_URL, params=params, timeout=30)
    resp.raise_for_status()

    papers, total = _parse_feed(resp.text)

    return SearchResult(
        papers=papers, total_count=total, start=start, max_results=max_results
    )


def get_papers(*arxiv_ids: str) -> list[Paper]:
    """Fetch specific papers by arXiv ID.

    Parameters
    ----------
    *arxiv_ids : str
        One or more arXiv IDs (e.g. "2103.12345", "hep-ex/0307015").

    Returns
    -------
    list[Paper]
    """
    if not arxiv_ids:
        raise ValueError("At least one arXiv ID is required")

    params = {
        "id_list": ",".join(arxiv_ids),
        "max_results": len(arxiv_ids),
    }

    _throttle()
    resp = requests.get(BASE_URL, params=params, timeout=30)
    resp.raise_for_status()

    papers, _ = _parse_feed(resp.text)
    return papers
