"""J-STAGE WebAPI client.

Provides typed access to the J-STAGE search API run by the Japan
Science and Technology Agency (JST). J-STAGE hosts ~3000 Japanese
scientific journals — including the flagship reactor physics
journals:

- **Journal of Nuclear Science and Technology (JNST)**
  ISSN 0022-3131 (print) / 1881-1248 (online), Atomic Energy
  Society of Japan. Primary outlet for Japanese cross-section
  evaluations (JENDL), MVP/SRAC benchmarks, BWR physics, fast-
  reactor work. Open access for most articles after a one-year
  embargo.
- **Transactions of the AESJ** (廃止 — old name).
- **Mechanical Engineering Journal** (Bridge of nuclear thermal-
  hydraulics work).
- **Journal of Power and Energy Systems**.

Also covers Japanese-only content rarely indexed in Scopus or
OpenAlex (early Yamamoto, Nakagawa, Chiba papers; AESJ proceedings
retroactively digitised).

The API returns Atom XML only — no JSON option exists.

No authentication required.

Usage::

    from tools.research.jstage import search, search_simple, JNST_ISSN

    # Search by free text + author + year range
    results = search_simple(
        text="MOC pin power",
        author="Yamamoto",
        year_from=2010,
    )

    # Restrict to JNST
    results = search_simple(
        text="resonance self-shielding",
        issn=JNST_ISSN,
    )

    # Raw search (any combination of fields)
    results = search(
        title="discrete ordinates",
        author="Chiba",
        keyword="cross section",
    )

    # Each Article has full bibliographic metadata
    for art in results.articles:
        print(art.summary())
        print("  DOI:", art.doi, "| URL:", art.url)
"""

from __future__ import annotations

import re
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import Optional

import requests

BASE_URL = "https://api.jstage.jst.go.jp/searchapi/do"

# Service codes (manual_api.pdf):
#   1 — Published volumes for a journal
#   2 — Article table-of-contents for a volume/issue
#   3 — Article search (keyword across articles) — what we use
#   4 — Journal search (find journals by title/ISSN/keyword)
SERVICE_ARTICLE_SEARCH = 3
SERVICE_JOURNAL_SEARCH = 4

# Common journals (ISSN — print). Use these with `issn=` for filtering.
JNST_ISSN = "0022-3131"  # Journal of Nuclear Science and Technology
JNST_ISSN_ONLINE = "1881-1248"
MEJ_ISSN = "2187-9745"  # Mechanical Engineering Journal
JPES_ISSN = "1881-3062"  # Journal of Power and Energy Systems

# J-STAGE uses a custom XML schema, NOT pure Atom. Most elements are
# unnamespaced; only PRISM bibliographic fields are namespaced. Each
# bilingual field has <en> and <ja> children — we prefer English.
_PRISM = "http://prismstandard.org/namespaces/basic/2.0/"
_NS = {
    "prism": _PRISM,
    "opensearch": "http://a9.com/-/spec/opensearch/1.1/",
}

_MIN_REQUEST_INTERVAL = 0.5
_last_request_time = 0.0


def _throttle() -> None:
    global _last_request_time
    elapsed = time.monotonic() - _last_request_time
    if elapsed < _MIN_REQUEST_INTERVAL:
        time.sleep(_MIN_REQUEST_INTERVAL - elapsed)
    _last_request_time = time.monotonic()


def _clean(text: Optional[str]) -> str:
    if text is None:
        return ""
    return re.sub(r"\s+", " ", text).strip()


def _local(tag: str) -> str:
    """Return the local part of a (possibly namespaced) tag."""
    return tag.rsplit("}", 1)[-1] if "}" in tag else tag


def _bilingual_text(elem: Optional[ET.Element], prefer: str = "en") -> str:
    """Pull text out of a J-STAGE bilingual <en>/<ja> wrapper.

    J-STAGE wraps human-readable strings as
    ``<field><en><![CDATA[..]]></en><ja>..</ja></field>``. Some
    fields skip the wrapper and put text directly in ``<field>``.
    """
    if elem is None:
        return ""
    # Try preferred language first, then the other.
    order = (prefer, "ja" if prefer == "en" else "en")
    for lang in order:
        for child in elem:
            if _local(child.tag) == lang:
                # Some <en> wrappers contain a <name> grandchild
                # (authors); recurse one level if no direct text.
                if child.text and child.text.strip():
                    return _clean(child.text)
                for gc in child:
                    if gc.text and gc.text.strip():
                        return _clean(gc.text)
    if elem.text and elem.text.strip():
        return _clean(elem.text)
    return ""


def _find_local(parent: ET.Element, name: str) -> Optional[ET.Element]:
    """First child matching ``name`` (local-name match, ignores ns)."""
    for child in parent:
        if _local(child.tag) == name:
            return child
    return None


def _find_all_local(parent: ET.Element, name: str) -> list[ET.Element]:
    return [c for c in parent if _local(c.tag) == name]


# ── Data classes ──────────────────────────────────────────────────────


@dataclass
class Article:
    """A single J-STAGE article record."""

    title: str = ""
    authors: list[str] = field(default_factory=list)
    journal_title: str = ""
    issn: str = ""
    eissn: str = ""
    volume: str = ""
    issue: str = ""
    start_page: str = ""
    end_page: str = ""
    publication_year: str = ""
    publication_date: str = ""
    doi: str = ""
    abstract: str = ""
    keywords: list[str] = field(default_factory=list)
    url: str = ""  # J-STAGE landing page
    pdf_url: str = ""
    material: str = ""  # journal code (cdjournal)

    @classmethod
    def from_entry(cls, entry: ET.Element) -> Article:
        """Parse from a J-STAGE <entry> element.

        J-STAGE's "Atom-like" feed actually carries non-namespaced
        custom elements with ``<en>`` / ``<ja>`` bilingual wrappers
        plus a few PRISM-namespaced bibliographic fields. This parser
        prefers English text and falls back to Japanese when only the
        Japanese form is populated.
        """
        title = _bilingual_text(_find_local(entry, "article_title"))
        # Some feeds also expose a flat <title> element — use as backup.
        if not title:
            t = _find_local(entry, "title")
            if t is not None and t.text:
                title = _clean(t.text)

        # Authors: <author><en><name>...</name></en><ja>...</ja></author>
        authors: list[str] = []
        for a in _find_all_local(entry, "author"):
            for lang_elem in a:
                if _local(lang_elem.tag) not in ("en", "ja"):
                    continue
                for sub in lang_elem:
                    if _local(sub.tag) == "name" and sub.text and sub.text.strip():
                        authors.append(_clean(sub.text))
                        break
                else:
                    if lang_elem.text and lang_elem.text.strip():
                        authors.append(_clean(lang_elem.text))
                if authors and authors[-1]:  # one author per <author>
                    break

        # Links: <article_link><en>URL</en></article_link>
        url = _bilingual_text(_find_local(entry, "article_link"))
        pdf = ""
        # Sometimes a <pdf_link> sibling exists.
        pl = _find_local(entry, "article_pdflink") or _find_local(entry, "pdf_link")
        if pl is not None:
            pdf = _bilingual_text(pl)

        # PRISM-namespaced fields are reliable.
        def _prism(name: str) -> str:
            t = entry.findtext(f"prism:{name}", default="", namespaces=_NS)
            return _clean(t) if t else ""

        # Plain unnamespaced bibliographic fields.
        def _plain(name: str) -> str:
            el = _find_local(entry, name)
            return _clean(el.text) if el is not None and el.text else ""

        journal = _bilingual_text(_find_local(entry, "material_title")) or _prism("publicationName")
        issn = _prism("issn")
        eissn = _prism("eIssn")
        volume = _prism("volume")
        issue = _prism("number")
        start_page = _prism("startingPage")
        end_page = _prism("endingPage")
        pubyear = _plain("pubyear")
        pubdate = _prism("publicationDate") or _plain("published") or _plain("updated")
        doi = _prism("doi") or _plain("doi")
        abstract = _bilingual_text(_find_local(entry, "article_abstract")) or _plain("abstract")
        material = _plain("cdjournal")

        kw_elem = _find_local(entry, "article_keyword") or _find_local(entry, "keyword")
        kw_raw = _bilingual_text(kw_elem) if kw_elem is not None else ""
        keywords = [k.strip() for k in re.split(r"[,;／、]", kw_raw) if k.strip()] if kw_raw else []

        return cls(
            title=title,
            authors=authors,
            journal_title=journal,
            issn=issn,
            eissn=eissn,
            volume=volume,
            issue=issue,
            start_page=start_page,
            end_page=end_page,
            publication_year=pubyear or (pubdate[:4] if pubdate else ""),
            publication_date=pubdate,
            doi=doi.replace("https://doi.org/", "").replace("http://doi.org/", ""),
            abstract=abstract,
            keywords=keywords,
            url=url,
            pdf_url=pdf,
            material=material,
        )

    @property
    def first_author(self) -> str:
        return self.authors[0] if self.authors else "Unknown"

    def summary(self) -> str:
        """One-line summary for display."""
        year = self.publication_year or "n.d."
        loc = ""
        if self.journal_title:
            loc = f" — {self.journal_title}"
            if self.volume:
                loc += f" {self.volume}"
                if self.issue:
                    loc += f"({self.issue})"
            if self.start_page:
                loc += f":{self.start_page}"
        return f"[J-STAGE] {self.first_author} ({year}): {self.title}{loc}"


@dataclass
class SearchResult:
    """Container for paginated J-STAGE search results."""

    articles: list[Article]
    total_count: int
    start: int
    count: int

    @property
    def total_pages(self) -> int:
        if self.count <= 0:
            return 0
        return (self.total_count + self.count - 1) // self.count


# ── XML parser ────────────────────────────────────────────────────────


class JStageAPIError(RuntimeError):
    """Raised when J-STAGE returns a non-zero <status> code."""


def _parse_feed(xml_text: str) -> tuple[list[Article], int]:
    """Parse a J-STAGE search feed into Articles and total result count.

    The feed uses ``opensearch:totalResults`` for the count and a
    flat list of ``<entry>`` children for the records (these
    ``entry`` elements are unnamespaced — not Atom).

    Raises
    ------
    JStageAPIError
        If the feed's ``<result><status>`` element reports an
        ``ERR_xxx`` code (the most common cause is an invalid
        parameter combination — see :func:`search`).
    """
    root = ET.fromstring(xml_text)

    # Surface API-level errors as exceptions.
    result = next((c for c in root if _local(c.tag) == "result"), None)
    if result is not None:
        status = next((c for c in result if _local(c.tag) == "status"), None)
        if status is not None and status.text and status.text.startswith("ERR"):
            msg = next((c for c in result if _local(c.tag) == "message"), None)
            raise JStageAPIError(
                f"J-STAGE returned {status.text}"
                + (f": {msg.text}" if msg is not None and msg.text else "")
            )

    total_str = root.findtext(
        "opensearch:totalResults", default="", namespaces=_NS
    )
    if not total_str:
        for child in root.iter():
            if _local(child.tag) == "totalResults" and child.text:
                total_str = child.text
                break
    try:
        total = int(re.sub(r"\D", "", total_str)) if total_str else 0
    except ValueError:
        total = 0

    entries = [c for c in root.iter() if _local(c.tag) == "entry"]
    articles = [Article.from_entry(e) for e in entries]
    return articles, total


# ── API functions ─────────────────────────────────────────────────────


def search(
    *,
    text: Optional[str] = None,
    title: Optional[str] = None,
    author: Optional[str] = None,
    affil: Optional[str] = None,
    keyword: Optional[str] = None,
    abstract: Optional[str] = None,
    issn: Optional[str] = None,
    cdjournal: Optional[str] = None,
    pubyearfrom: Optional[int] = None,
    pubyearto: Optional[int] = None,
    sortflg: int = 1,
    count: int = 20,
    start: int = 1,
    lang: str = "en",
) -> SearchResult:
    """Search J-STAGE articles.

    All search-term parameters AND-combine. Pass ``None`` (or omit)
    to skip a clause. At least one search term must be provided.

    Parameters
    ----------
    text : str
        Free-text search across title, abstract, keyword, and
        author. Most useful starting point.
    title : str
        Match within article title.
    author : str
        Match within author name. Western names work — prefer the
        romaji form (e.g. ``"Yamamoto"``, ``"Chiba"``).
    affil : str
        Author affiliation (institution).
    keyword : str
        Match within article keywords.
    abstract : str
        Match within abstract text.
    issn : str
        Restrict to a single journal by ISSN. Convenience constants:
        :data:`JNST_ISSN`, :data:`JNST_ISSN_ONLINE`, :data:`MEJ_ISSN`,
        :data:`JPES_ISSN`.
    cdjournal : str
        J-STAGE journal code (alternative to ISSN). Visible in the
        URL path when browsing a journal on j-stage.jst.go.jp.
    pubyearfrom, pubyearto : int
        Inclusive year bounds. Either or both.
    sortflg : int
        1 = sort by relevance (default); 2 = newest publication date.
    count : int
        Results per page (max 1000).
    start : int
        1-indexed result offset (NOT zero-indexed — J-STAGE is unusual).
    lang : str
        Response language hint: ``"ja"`` or ``"en"``.

    Returns
    -------
    SearchResult

    Raises
    ------
    JStageAPIError
        If J-STAGE rejects the parameter combination (status
        ``ERR_001``) or an individual parameter is malformed (e.g.
        ``ERR_008`` for an unhyphenated ISSN).

    Notes
    -----
    J-STAGE's article search is finicky about which parameters can
    appear together. The empirical rule: you may either use
    ``text`` (broad full-text) plus **one** narrowing field, OR
    combine **two** specific fields without ``text``. Mixing
    ``text`` with two specific fields raises ``ERR_001``.

    Validated 2-parameter combinations:

    ✅ ``text + author``
    ✅ ``text + issn``
    ✅ ``text + pubyearfrom`` (and/or ``pubyearto``)
    ✅ ``author + issn``
    ✅ ``title + author``

    Known to fail:

    ❌ ``text + author + issn`` — returns ``ERR_001``
    ❌ ``keyword + issn`` — returns ``ERR_001``
    ❌ ``issn + pubyearfrom`` — returns ``ERR_001``
    ❌ ``cdjournal`` alone — returns ``ERR_001``

    When in doubt, search broadly with one field-specific filter
    and post-filter the returned :class:`Article` list locally on
    ``issn`` / ``pubyear`` / ``journal_title``. The ISSN must
    include the hyphen (``"0022-3131"``, not ``"00223131"``).
    """
    params: dict = {"service": SERVICE_ARTICLE_SEARCH, "lang": lang}
    if text is not None:
        params["text"] = text
    if title is not None:
        params["title"] = title
    if author is not None:
        params["author"] = author
    if affil is not None:
        params["affil"] = affil
    if keyword is not None:
        params["keyword"] = keyword
    if abstract is not None:
        params["abstract"] = abstract
    if issn is not None:
        params["issn"] = issn
    if cdjournal is not None:
        params["cdjournal"] = cdjournal
    if pubyearfrom is not None:
        params["pubyearfrom"] = pubyearfrom
    if pubyearto is not None:
        params["pubyearto"] = pubyearto
    params["sortflg"] = sortflg
    params["count"] = count
    params["start"] = start

    has_term = any(
        params.get(k) for k in ("text", "title", "author", "affil",
                                "keyword", "abstract", "issn", "cdjournal")
    )
    if not has_term and pubyearfrom is None and pubyearto is None:
        raise ValueError("Provide at least one search term or year filter.")

    _throttle()
    resp = requests.get(BASE_URL, params=params, timeout=30)
    resp.raise_for_status()

    articles, total = _parse_feed(resp.text)
    return SearchResult(
        articles=articles, total_count=total, start=start, count=count
    )


# Alias kept for parity with the other research modules. ``search``
# already accepts kwargs cleanly, but ``search_simple`` makes the
# intent explicit when called via ``-c "..."``.
search_simple = search


def search_journal(
    *,
    text: Optional[str] = None,
    title: Optional[str] = None,
    issn: Optional[str] = None,
    cdjournal: Optional[str] = None,
    count: int = 20,
    start: int = 1,
    lang: str = "en",
) -> str:
    """Search the J-STAGE journal directory (service=4).

    Use this to discover ``cdjournal`` codes or ISSN values when you
    only know a journal's name. Returns the raw Atom XML — for the
    handful of fields you typically need it's easier than building a
    full ``Journal`` dataclass.

    At least one of ``text``, ``title``, ``issn``, or ``cdjournal``
    must be provided.
    """
    params: dict = {"service": SERVICE_JOURNAL_SEARCH, "lang": lang,
                    "count": count, "start": start}
    if text is not None:
        params["text"] = text
    if title is not None:
        params["title"] = title
    if issn is not None:
        params["issn"] = issn
    if cdjournal is not None:
        params["cdjournal"] = cdjournal
    if not any(params.get(k) for k in ("text", "title", "issn", "cdjournal")):
        raise ValueError("Provide at least one journal search term.")

    _throttle()
    resp = requests.get(BASE_URL, params=params, timeout=30)
    resp.raise_for_status()
    return resp.text
