"""IAEA INIS (International Nuclear Information System) API client.

Provides typed access to the INIS InvenioRDM REST API for searching
the world's largest nuclear-specific literature database (4M+ records).

Covers journal articles, technical reports, conference papers, theses,
and standards across all nuclear topics. Indexes publications not found
in CrossRef or OpenAlex (IAEA reports, Eastern European journals, etc.).

No authentication required. Returns metadata only (not full text).

Usage::

    from tools.research.inis import search, get_record

    # Free-text search
    results = search("thermal hydraulics LOCA")

    # Filtered search
    results = search("collision probability", size=5, sort="newest")

    # Fetch a specific record
    record = get_record("12345678")
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Optional

import requests

BASE_URL = "https://inis.iaea.org/api/records"

_MIN_REQUEST_INTERVAL = 0.5
_last_request_time = 0.0


def _throttle() -> None:
    global _last_request_time
    elapsed = time.monotonic() - _last_request_time
    if elapsed < _MIN_REQUEST_INTERVAL:
        time.sleep(_MIN_REQUEST_INTERVAL - elapsed)
    _last_request_time = time.monotonic()


# ── Data classes ──────────────────────────────────────────────────────


@dataclass
class Record:
    """A single INIS bibliographic record."""

    inis_id: str
    title: str
    authors: list[str] = field(default_factory=list)
    description: str = ""
    publication_date: str = ""
    resource_type: str = ""
    subjects: list[str] = field(default_factory=list)
    identifiers: dict[str, str] = field(default_factory=dict)
    publisher: str = ""
    language: str = ""
    url: str = ""

    @classmethod
    def from_json(cls, data: dict) -> Record:
        """Parse a record from the INIS API JSON response."""
        metadata = data.get("metadata", {})

        # Title
        title = metadata.get("title", "")

        # Authors — creators list
        authors = []
        for creator in metadata.get("creators", []):
            person = creator.get("person_or_org", {})
            name = person.get("name", "")
            if not name:
                given = person.get("given_name", "")
                family = person.get("family_name", "")
                name = f"{given} {family}".strip()
            if name:
                authors.append(name)

        # Description / abstract
        description = ""
        for desc in metadata.get("descriptions", []):
            if desc.get("type", {}).get("id") == "abstract" or not description:
                description = desc.get("description", "")

        # Publication date
        pub_date = metadata.get("publication_date", "")

        # Resource type
        resource_type = ""
        rt = metadata.get("resource_type", {})
        resource_type = rt.get("title", {}).get("en", rt.get("id", ""))

        # Subjects
        subjects = []
        for subj in metadata.get("subjects", []):
            s = subj.get("subject", "")
            if s:
                subjects.append(s)

        # Identifiers (DOI, INIS ref number, etc.)
        identifiers = {}
        for ident in metadata.get("identifiers", []):
            scheme = ident.get("scheme", "other")
            value = ident.get("identifier", "")
            if value:
                identifiers[scheme] = value

        # Publisher
        publisher = metadata.get("publisher", "")

        # Language
        languages = metadata.get("languages", [])
        language = languages[0].get("id", "") if languages else ""

        # Self link
        url = data.get("links", {}).get("self_html", "")

        return cls(
            inis_id=str(data.get("id", "")),
            title=title,
            authors=authors,
            description=description,
            publication_date=pub_date,
            resource_type=resource_type,
            subjects=subjects,
            identifiers=identifiers,
            publisher=publisher,
            language=language,
            url=url,
        )

    @property
    def doi(self) -> str:
        return self.identifiers.get("doi", "")

    def summary(self) -> str:
        """One-line summary for display."""
        year = self.publication_date[:4] if self.publication_date else "n.d."
        first_author = self.authors[0].split(",")[0] if self.authors else "Unknown"
        return f"[INIS:{self.inis_id}] {first_author} ({year}): {self.title}"


@dataclass
class SearchResult:
    """Container for paginated INIS search results."""

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
    q: str,
    *,
    size: int = 10,
    page: int = 1,
    sort: str = "newest",
) -> SearchResult:
    """Search INIS for nuclear literature.

    Parameters
    ----------
    q : str
        Search query string (supports Elasticsearch query syntax).
    size : int
        Results per page (default 10).
    page : int
        Page number (1-indexed, default 1).
    sort : str
        Sort order: "newest", "oldest", or "bestmatch" (default "newest").

    Returns
    -------
    SearchResult
    """
    params: dict = {
        "q": q,
        "size": size,
        "page": page,
        "sort": sort,
    }

    _throttle()
    resp = requests.get(BASE_URL, params=params, timeout=30)
    resp.raise_for_status()

    data = resp.json()
    hits = data.get("hits", {})
    total = hits.get("total", 0)
    if isinstance(total, dict):
        total = total.get("value", 0)

    records = [Record.from_json(item) for item in hits.get("hits", [])]

    return SearchResult(records=records, total_count=total, page=page, size=size)


def get_record(record_id: str) -> Record:
    """Fetch a single INIS record by ID.

    Parameters
    ----------
    record_id : str
        The INIS record identifier.

    Returns
    -------
    Record
    """
    _throttle()
    resp = requests.get(f"{BASE_URL}/{record_id}", timeout=30)
    resp.raise_for_status()
    return Record.from_json(resp.json())
