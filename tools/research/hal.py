"""HAL Science (archives-ouvertes.fr) API client.

Provides typed access to the HAL Solr-based search API, the French
national open archive for scholarly publications. Covers journal
articles, conference papers, theses (HDR + PhD), reports, and book
chapters across all disciplines, with very strong coverage of
French nuclear-engineering output (CEA, IRSN, EDF R&D, IRFU, CNRS
LPSC/LP2I/IJCLab, Subatech).

Useful complement to OSTI / INIS / OpenAlex when looking for:

- CEA reports and theses (often deposited in HAL but absent from
  OSTI and behind paywalls in journal-only sources).
- French-language nuclear engineering and reactor physics work
  (e.g. Sanchez, Hébert, Reuss, Coste-Delclaux).
- Open-access PDFs of conference papers (PHYSOR, M&C, ICAPP,
  ICONE, NURETH) when authors deposit accepted manuscripts.
- HDR theses ("Habilitation à diriger des recherches") that act
  as comprehensive monographs on a researcher's career and rarely
  appear in other databases.

No authentication required. Default polite throttle 0.5 s.

Usage::

    from tools.research.hal import search, get_record, get_bibtex

    # Free-text search
    results = search("collision probability cylindrical")

    # Fielded search using helper kwargs
    results = search_simple(
        title="probabilités de collision",
        author="Sanchez",
        year_from=1980,
    )

    # Native Solr query (full power)
    results = search('authFullName_s:"Richard Sanchez" AND docType_s:THESE')

    # Filter by document type (journal article, conference, thesis, ...)
    results = search("Hébert", doc_type="ART", rows=20)

    # Resolve a known record by HAL id or DOI
    rec = get_record("hal-01242305")
    rec = get_record(doi="10.13182/NSE08-64")

    # BibTeX export
    bib = get_bibtex("hal-01242305")
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Optional, Union

import requests

BASE_URL = "https://api.archives-ouvertes.fr/search"

# Default field selection — keeps payloads small while exposing
# everything Record.from_json knows how to read.
_DEFAULT_FIELDS = (
    "halId_s,docid,doiId_s,title_s,subTitle_s,abstract_s,"
    "authFullName_s,authLastName_s,authIdHal_s,"
    "journalTitle_s,journalIssn_s,volume_s,issue_s,page_s,"
    "producedDate_tdate,producedDateY_i,publicationDate_tdate,"
    "docType_s,language_s,domainAllCode_s,keyword_s,"
    "uri_s,fileMain_s,files_s,openAccess_bool,"
    "conferenceTitle_s,city_s,country_s,publisher_s"
)

# docType_s vocabulary (most common values; full list at
# https://api.archives-ouvertes.fr/ref/doctype/).
DOC_TYPES = {
    "ART": "Journal article",
    "COMM": "Conference paper",
    "POSTER": "Conference poster",
    "PROCEEDINGS": "Conference proceedings",
    "OUV": "Book",
    "COUV": "Book chapter",
    "DOUV": "Edited volume / proceedings book",
    "REPORT": "Report (technical or research)",
    "THESE": "PhD thesis",
    "HDR": "Habilitation thesis",
    "MEM": "Master's thesis",
    "PATENT": "Patent",
    "PRESCONF": "Conference invited talk",
    "SOFTWARE": "Software",
    "PREPRINT": "Preprint",
    "OTHER": "Other",
    "UNDEFINED": "Unspecified",
}

_MIN_REQUEST_INTERVAL = 0.5
_last_request_time = 0.0


def _throttle() -> None:
    global _last_request_time
    elapsed = time.monotonic() - _last_request_time
    if elapsed < _MIN_REQUEST_INTERVAL:
        time.sleep(_MIN_REQUEST_INTERVAL - elapsed)
    _last_request_time = time.monotonic()


def _first(value):
    """Solr returns multi-valued fields as lists even for single hits."""
    if isinstance(value, list):
        return value[0] if value else ""
    return value or ""


def _as_list(value) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [v for v in value if v]
    return [value] if value else []


# ── Data classes ──────────────────────────────────────────────────────


@dataclass
class Record:
    """A single HAL bibliographic record."""

    hal_id: str
    docid: str = ""
    title: str = ""
    subtitle: str = ""
    authors: list[str] = field(default_factory=list)
    author_last_names: list[str] = field(default_factory=list)
    author_hal_ids: list[str] = field(default_factory=list)
    abstract: str = ""
    doi: str = ""
    journal: str = ""
    journal_issn: str = ""
    volume: str = ""
    issue: str = ""
    pages: str = ""
    publication_date: str = ""  # "produced" date — when the work appeared
    publication_year: int = 0
    deposit_date: str = ""  # when uploaded to HAL
    doc_type: str = ""
    language: str = ""
    domains: list[str] = field(default_factory=list)
    keywords: list[str] = field(default_factory=list)
    url: str = ""  # HAL landing page
    pdf_url: str = ""  # main full-text file
    files: list[str] = field(default_factory=list)  # all attached files
    is_open_access: bool = False
    conference_title: str = ""
    conference_city: str = ""
    conference_country: str = ""
    publisher: str = ""

    @classmethod
    def from_json(cls, data: dict) -> Record:
        """Parse a record from a HAL Solr-style search hit."""
        title = _first(data.get("title_s"))
        subtitle = _first(data.get("subTitle_s"))
        abstract = _first(data.get("abstract_s"))

        # Date fields — HAL may return either ISO timestamps or the
        # year alone depending on what the depositor entered.
        produced = _first(data.get("producedDate_tdate"))
        deposit = _first(data.get("publicationDate_tdate"))
        year = data.get("producedDateY_i")
        if isinstance(year, list):
            year = year[0] if year else 0
        try:
            year = int(year) if year else 0
        except (TypeError, ValueError):
            year = 0

        return cls(
            hal_id=_first(data.get("halId_s")),
            docid=str(data.get("docid", "") or ""),
            title=title,
            subtitle=subtitle,
            authors=_as_list(data.get("authFullName_s")),
            author_last_names=_as_list(data.get("authLastName_s")),
            author_hal_ids=_as_list(data.get("authIdHal_s")),
            abstract=abstract,
            doi=_first(data.get("doiId_s")),
            journal=_first(data.get("journalTitle_s")),
            journal_issn=_first(data.get("journalIssn_s")),
            volume=_first(data.get("volume_s")),
            issue=_first(data.get("issue_s")),
            pages=_first(data.get("page_s")),
            publication_date=produced,
            publication_year=year,
            deposit_date=deposit,
            doc_type=_first(data.get("docType_s")),
            language=_first(data.get("language_s")),
            domains=_as_list(data.get("domainAllCode_s")),
            keywords=_as_list(data.get("keyword_s")),
            url=_first(data.get("uri_s")),
            pdf_url=_first(data.get("fileMain_s")),
            files=_as_list(data.get("files_s")),
            is_open_access=bool(_first(data.get("openAccess_bool")) or False),
            conference_title=_first(data.get("conferenceTitle_s")),
            conference_city=_first(data.get("city_s")),
            conference_country=_first(data.get("country_s")),
            publisher=_first(data.get("publisher_s")),
        )

    @property
    def first_author(self) -> str:
        return self.authors[0] if self.authors else "Unknown"

    @property
    def doc_type_label(self) -> str:
        return DOC_TYPES.get(self.doc_type, self.doc_type or "Unknown")

    def summary(self) -> str:
        """One-line summary for display."""
        year = self.publication_year or (
            self.publication_date[:4] if self.publication_date else "n.d."
        )
        oa = " [OA]" if self.is_open_access or self.pdf_url else ""
        kind = f" [{self.doc_type}]" if self.doc_type else ""
        return f"[HAL:{self.hal_id}] {self.first_author} ({year}){kind}: {self.title}{oa}"


@dataclass
class SearchResult:
    """Container for paginated HAL search results."""

    records: list[Record]
    total_count: int
    start: int
    rows: int

    @property
    def page(self) -> int:
        if self.rows <= 0:
            return 1
        return self.start // self.rows + 1

    @property
    def total_pages(self) -> int:
        if self.rows <= 0:
            return 0
        return (self.total_count + self.rows - 1) // self.rows


# ── API functions ─────────────────────────────────────────────────────


def search(
    q: str,
    *,
    doc_type: Optional[Union[str, list[str]]] = None,
    year_from: Optional[int] = None,
    year_to: Optional[int] = None,
    open_access_only: bool = False,
    fl: str = _DEFAULT_FIELDS,
    rows: int = 10,
    start: int = 0,
    sort: str = "producedDate_tdate desc",
    portal: Optional[str] = None,
) -> SearchResult:
    """Search HAL with a Solr query.

    Parameters
    ----------
    q : str
        Solr query. Plain text is a free-text search; field-prefixed
        terms (e.g. ``authFullName_s:Sanchez``) target a specific
        Solr field. Use parentheses and quoting for multi-word values.
    doc_type : str or list[str], optional
        Restrict to one or more docType_s codes (see
        :data:`DOC_TYPES`). Passes through as an ``fq`` filter.
    year_from, year_to : int, optional
        Inclusive bounds on ``producedDateY_i``. Either or both may
        be given.
    open_access_only : bool
        If True, restrict to records flagged ``openAccess_bool:true``.
    fl : str
        Comma-separated Solr field list. Defaults cover everything
        :class:`Record` parses; widen only if you need extras.
    rows : int
        Page size (max 10000 per HAL docs). Default 10.
    start : int
        Zero-indexed offset into the result set. Default 0.
    sort : str
        Solr sort clause. Default sorts newest-first by produced date.
        Other useful values: ``"score desc"`` (relevance),
        ``"submittedDate_tdate desc"``.
    portal : str, optional
        Restrict the query to a HAL portal (lowercase, e.g. ``"cea"``,
        ``"irsn"``, ``"in2p3"``) or collection (uppercase code).

    Returns
    -------
    SearchResult

    Notes
    -----
    Solr field reference (selected — full list under HAL ref docs):

    ============================  =====================================
    Field                         Meaning
    ============================  =====================================
    ``title_s``                   Title (string, exact match)
    ``title_t``                   Title (analyzed text — better recall)
    ``authFullName_s``            "First Last" exact-match string
    ``authFullName_t``            "First Last" analyzed (diacritic-folded)
    ``authLastName_s``            Surname (exact match)
    ``authLastName_t``            Surname (analyzed, diacritic-folded)
    ``authIdHal_s``               HAL author profile id
    ``producedDateY_i``           Year (int) of publication
    ``producedDate_tdate``        Full ISO date of publication
    ``docType_s``                 Document type code (see DOC_TYPES)
    ``journalTitle_s``            Journal title (exact match)
    ``journalTitle_t``            Journal title (analyzed)
    ``journalIssn_s``             Journal ISSN
    ``keyword_s``                 Author-supplied keyword
    ``abstract_s``                Abstract text
    ``doiId_s``                   DOI (lowercased)
    ``halId_s``                   HAL deposit id (e.g. "hal-01242305")
    ``domainAllCode_s``           Discipline code (e.g. "phys.nucl-ex")
    ``openAccess_bool``           True iff a full-text file is attached
    ============================  =====================================

    The ``_t`` (text) variants are tokenized and diacritic-folded —
    use them whenever the user-supplied query may differ in case or
    accent from the deposited record (the common case for French
    surnames). Use ``_s`` (string) for exact-match filters such as
    HAL ids, DOIs, or controlled vocabularies.
    """
    endpoint = f"{BASE_URL}/{portal}/" if portal else f"{BASE_URL}/"

    params: dict = {
        "q": q,
        "wt": "json",
        "fl": fl,
        "rows": rows,
        "start": start,
        "sort": sort,
    }

    fq: list[str] = []
    if doc_type is not None:
        types = [doc_type] if isinstance(doc_type, str) else list(doc_type)
        if len(types) == 1:
            fq.append(f"docType_s:{types[0]}")
        else:
            fq.append("docType_s:(" + " OR ".join(types) + ")")
    if year_from is not None or year_to is not None:
        lo = year_from if year_from is not None else "*"
        hi = year_to if year_to is not None else "*"
        fq.append(f"producedDateY_i:[{lo} TO {hi}]")
    if open_access_only:
        fq.append("openAccess_bool:true")
    if fq:
        params["fq"] = fq  # requests will repeat the param for list values

    _throttle()
    resp = requests.get(endpoint, params=params, timeout=30)
    resp.raise_for_status()

    data = resp.json()
    response = data.get("response", {})
    total = response.get("numFound", 0)
    docs = response.get("docs", [])

    records = [Record.from_json(item) for item in docs]
    return SearchResult(records=records, total_count=total, start=start, rows=rows)


def search_simple(
    *,
    keywords: Optional[str] = None,
    title: Optional[str] = None,
    author: Optional[str] = None,
    journal: Optional[str] = None,
    year_from: Optional[int] = None,
    year_to: Optional[int] = None,
    doc_type: Optional[Union[str, list[str]]] = None,
    open_access_only: bool = False,
    rows: int = 10,
    start: int = 0,
    sort: str = "producedDate_tdate desc",
    portal: Optional[str] = None,
) -> SearchResult:
    """Convenience wrapper that builds a Solr ``q`` from kwargs.

    All kwargs are AND-combined. Multi-word values are quoted. Pass
    ``None`` (or omit) to skip a clause.

    Examples
    --------
    >>> search_simple(author="Sanchez", title="probabilités collision")
    >>> search_simple(keywords="MOC", journal="Annals of Nuclear Energy",
    ...               year_from=2010, doc_type="ART")
    """
    clauses: list[str] = []

    def _quote(value: str) -> str:
        v = value.strip()
        if " " in v or any(c in v for c in ":()[]{}\"'"):
            v = v.replace('"', r"\"")
            return f'"{v}"'
        return v

    if keywords:
        clauses.append(_quote(keywords))
    if title:
        clauses.append(f"title_t:{_quote(title)}")
    if author:
        # Use the analyzed-text fields: tokenized + diacritic-folded,
        # so "Hebert" matches "Hébert" and "Alain Hebert" matches
        # the full "Alain Hébert" string regardless of token order.
        a = author.strip()
        field_name = "authFullName_t" if " " in a else "authLastName_t"
        clauses.append(f"{field_name}:{_quote(a)}")
    if journal:
        clauses.append(f"journalTitle_t:{_quote(journal)}")

    q = " AND ".join(clauses) if clauses else "*:*"

    return search(
        q,
        doc_type=doc_type,
        year_from=year_from,
        year_to=year_to,
        open_access_only=open_access_only,
        rows=rows,
        start=start,
        sort=sort,
        portal=portal,
    )


def get_record(
    hal_id: Optional[str] = None,
    *,
    doi: Optional[str] = None,
    docid: Optional[str] = None,
) -> Record:
    """Fetch a single HAL record by HAL id, DOI, or internal docid.

    Exactly one of the identifiers must be provided.

    Parameters
    ----------
    hal_id : str, optional
        HAL deposit id, e.g. ``"hal-01242305"`` or ``"tel-00012345"``.
    doi : str, optional
        DOI (with or without the ``https://doi.org/`` prefix).
    docid : str, optional
        Internal HAL document id (numeric).

    Returns
    -------
    Record

    Raises
    ------
    LookupError
        If the lookup returns no documents.
    ValueError
        If no identifier (or more than one) is supplied.
    """
    given = [x for x in (hal_id, doi, docid) if x]
    if len(given) != 1:
        raise ValueError("Provide exactly one of hal_id, doi, or docid.")

    if hal_id:
        q = f"halId_s:{hal_id}"
    elif docid:
        q = f"docid:{docid}"
    else:
        d = doi.replace("https://doi.org/", "").replace("http://doi.org/", "")
        # DOIs may contain Solr-special characters — escape them.
        for ch in r":/()[]{}":
            d = d.replace(ch, f"\\{ch}")
        q = f"doiId_s:{d}"

    result = search(q, rows=1, sort="score desc")
    if not result.records:
        raise LookupError(f"No HAL record matched query: {q}")
    return result.records[0]


def get_bibtex(hal_id: str) -> str:
    """Fetch a BibTeX entry for a HAL record.

    Uses the search API with ``wt=bibtex``, which is the documented
    structured-export route. The web landing page at
    ``hal.science/<hal_id>/bibtex`` returns the HTML page, not the
    raw entry.

    Parameters
    ----------
    hal_id : str
        HAL deposit id, e.g. ``"hal-01242305"``.

    Returns
    -------
    str
        BibTeX entry text.
    """
    _throttle()
    resp = requests.get(
        f"{BASE_URL}/",
        params={"q": f"halId_s:{hal_id}", "wt": "bibtex", "rows": 1},
        timeout=30,
    )
    resp.raise_for_status()
    return resp.text
