---
name: research
description: >
  Use when the user or an agent needs to search scientific literature,
  resolve DOIs, explore citation graphs, or query nuclear data.
  8 databases: OSTI, arXiv, Scopus, INIS, OpenAlex, CrossRef, Semantic
  Scholar, IAEA-NDS. Examples: "Find papers on collision probability",
  "Look up Bailey 2009", "What cites this DOI?", "U-235 ground state data"
---

# Literature & Nuclear Data Search Skill

Search scientific databases for nuclear engineering literature, resolve
DOIs, explore citation networks, and query nuclear structure data.

## Available Databases

### Literature Search

| Database | Best for | Auth | Module |
|----------|----------|------|--------|
| **OSTI.gov** | DOE technical reports, national lab pubs | None | `tools.research.osti` |
| **arXiv** | Preprints, open-access papers | None | `tools.research.arxiv` |
| **Scopus** | Journal articles, citation counts | API key | `tools.research.elsevier` |
| **INIS** | Nuclear-specific literature (4M+ records, IAEA reports) | None | `tools.research.inis` |
| **OpenAlex** | Scholarly graph, OA PDFs, citation networks | None | `tools.research.openalex` |
| **CrossRef** | DOI resolution, journal metadata | None | `tools.research.crossref` |
| **Semantic Scholar** | AI-enhanced search, influential citations | None (key optional) | `tools.research.semantic_scholar` |

### Nuclear Data

| Database | Best for | Auth | Module |
|----------|----------|------|--------|
| **IAEA-NDS LiveChart** | Nuclear structure: half-lives, binding energies, gammas, levels | None | `tools.research.iaea_nds` |

## When to Use Which

- **DOE lab report or NUREG** → OSTI
- **Recent/unpublished work** → arXiv
- **Specific journal article** (by DOI, author, journal) → Scopus or CrossRef
- **IAEA reports, Eastern European journals** → INIS
- **Citation graph** ("what cites X?" / "what does X cite?") → OpenAlex or Semantic Scholar
- **Open-access PDF link** → OpenAlex (`oa_url` field)
- **Impact / influential citations** → Semantic Scholar (`influential_citation_count`)
- **DOI metadata resolution** → CrossRef (`get_work`)
- **Journal-specific search** → CrossRef (`journal_search` by ISSN)
- **Broad topic search** → search OSTI + OpenAlex + INIS in parallel
- **Nuclear structure data** (half-life, binding energy, gammas) → IAEA-NDS
- **BibTeX** → OSTI (`get_bibtex`) or arXiv (`bibtex_key`)

## Quick Reference

### OSTI.gov

```python
from tools.research.osti import search, get_record, get_bibtex

results = search("discrete ordinates cylindrical geometry")
results = search(title="collision probability", author="Hébert", rows=10)
record = get_record("1893820")
bibtex = get_bibtex("1893820")

record.title, record.authors, record.doi, record.fulltext_url, record.summary()
```

**Parameters**: `q`, `title`, `author`, `subject`, `fulltext`, `doi`, `publication_date_start`/`end` (MM/DD/YYYY), `has_fulltext`, `rows`, `page`, `sort`, `order`.

### arXiv

```python
from tools.research.arxiv import search, get_papers

results = search("neutron transport SN method")
results = search(title="discrete ordinates", author="Morel", category="nucl-th")
results = search(query="au:morel AND ti:discrete ordinates AND cat:nucl-th")
papers = get_papers("2103.12345", "hep-ex/0307015")

paper.title, paper.authors, paper.summary, paper.doi, paper.pdf_url
paper.journal_ref, paper.categories, paper.summary_line(), paper.bibtex_key()
```

**Field prefixes**: `ti:`, `au:`, `abs:`, `cat:`, `co:`, `jr:`, `all:`.
**Operators**: `AND`, `OR`, `ANDNOT`.
**Categories**: `nucl-th`, `nucl-ex`, `physics.comp-ph`, `math-ph`.

### Elsevier / Scopus

```python
from tools.research.elsevier import search, search_simple, get_abstract

results = search_simple(keywords="collision probability", author="Hébert",
                        journal="Annals of Nuclear Energy", year_from=2005)
results = search("TITLE-ABS-KEY(discrete ordinates) AND AUTH(morel) AND PUBYEAR > 2000")
results = search_simple(keywords="SN transport", sort="-citedby-count")
abstract = get_abstract(doi="10.1016/j.anucene.2023.01.001")

record.title, record.authors, record.doi, record.journal, record.cited_by
abstract.abstract, abstract.keywords, abstract.subject_areas
```

**Query codes**: `TITLE`, `ABS`, `KEY`, `TITLE-ABS-KEY`, `AUTH`, `DOI`, `SRCTITLE`, `PUBYEAR`, `DOCTYPE`, `OPENACCESS`.
**Requires API key** in `tools/research/elsevier.toml`.

### INIS (IAEA)

```python
from tools.research.inis import search, get_record

results = search("thermal hydraulics LOCA")
results = search("collision probability", size=5, sort="newest")
record = get_record("12345678")

record.title, record.authors, record.description, record.doi
record.resource_type, record.subjects, record.summary()
```

**Parameters**: `q` (Elasticsearch syntax), `size`, `page`, `sort` ("newest", "oldest", "bestmatch").

### OpenAlex

```python
from tools.research.openalex import search, get_work, get_citations, get_references

results = search("discrete ordinates neutron transport")
results = search("Monte Carlo", publication_year="2015-2023", open_access=True)
work = get_work("10.13182/NSE08-64")  # by DOI
citations = get_citations("10.13182/NSE08-64")  # who cites this?
references = get_references("10.13182/NSE08-64")  # what does it cite?

work.title, work.authors, work.abstract, work.doi, work.cited_by_count
work.is_open_access, work.oa_url, work.journal, work.concepts, work.summary()
```

**Filters**: `publication_year` (int or "2015-2023"), `type` ("article", "review", etc.), `open_access`.
**Sort**: "relevance_score:desc", "cited_by_count:desc", "publication_date:desc".

### CrossRef

```python
from tools.research.crossref import search, get_work, journal_search

work = get_work("10.13182/NSE08-64")  # DOI resolution
results = search("collision probability neutron")
results = search("SN transport", from_pub_date="2020", type="journal-article")
results = journal_search("0029-5639", query="discrete ordinates")  # by ISSN

work.title, work.authors, work.doi, work.container_title, work.published_date
work.is_referenced_by_count, work.abstract, work.summary()
```

**Filters**: `from_pub_date`, `until_pub_date`, `type`, `has_abstract`.
**Sort**: "relevance", "published", "is-referenced-by-count".

**Key ISSNs**:
```
NSE  = "0029-5639"   ANE  = "0306-4549"   NED  = "0029-5493"
NET  = "1738-5733"   JCP  = "0021-9991"   PNE  = "0149-1970"
JNM  = "0022-3115"   CPC  = "0010-4655"
```

### Semantic Scholar

```python
from tools.research.semantic_scholar import search, get_paper, get_citations, get_references

results = search("discrete ordinates neutron transport")
results = search("SN method", year="2015-2023", fields_of_study=["Physics"])
paper = get_paper("DOI:10.13182/NSE08-64")
paper = get_paper("ARXIV:2103.12345")
citations = get_citations("DOI:10.13182/NSE08-64")
references = get_references("DOI:10.13182/NSE08-64")

paper.title, paper.authors, paper.abstract, paper.doi, paper.year
paper.citation_count, paper.influential_citation_count
paper.is_open_access, paper.oa_pdf_url, paper.journal, paper.summary()
```

**ID formats**: `DOI:...`, `ARXIV:...`, `PMID:...`, `CorpusID:...`, or raw S2 hash.
**Filters**: `year` ("2020", "2015-2023", "2020-", "-2015"), `fields_of_study`, `open_access_only`.

### IAEA-NDS LiveChart (Nuclear Data)

```python
from tools.research.iaea_nds import ground_state, ground_states, ground_states_by_z
from tools.research.iaea_nds import gammas, levels

# Ground state properties
gs = ground_state("U235")
gs.half_life, gs.binding_energy, gs.spin_parity, gs.atomic_mass, gs.abundance

# Multiple nuclides
isotopes = ground_states(["U235", "U238", "Pu239"])

# All isotopes of an element
uranium = ground_states_by_z(92)

# Gamma transitions
co60_gammas = gammas("Co60")
for g in co60_gammas:
    g.energy, g.relative_intensity, g.multipolarity

# Nuclear energy levels
fe56_levels = levels("Fe56")
for l in fe56_levels:
    l.energy, l.spin_parity, l.half_life
```

**Note**: Fission yields are NOT available through this API. Use ENDF libraries.

## Workflows

### Comprehensive Literature Search

```
1. Identify terms and scope (topic, aspect, time range)
2. Search in parallel:
   - OSTI: DOE reports          - INIS: IAEA + global nuclear
   - arXiv: preprints           - OpenAlex: journals + OA PDFs
3. Deduplicate by DOI
4. Rank by relevance, flag OA availability and citation counts
```

### Find a Specific Paper

```
1. If DOI known → CrossRef get_work() for metadata
2. If author+year+journal → Scopus search_simple() or CrossRef journal_search()
3. If partial info → OpenAlex search() or Semantic Scholar search()
4. Return full citation + access links
```

### Citation Network Analysis

```python
# "What are the key papers in this area?"
from tools.research.openalex import search, get_citations

# Find seed paper
results = search("Bailey discrete ordinates curvilinear", per_page=1)
seed = results.works[0]

# Who cites it?
citing = get_citations(seed.doi, per_page=50)
for w in citing.works:
    print(w.summary())  # includes citation count
```

```python
# "What are the influential citations?"
from tools.research.semantic_scholar import get_paper, get_citations

paper = get_paper("DOI:10.13182/NSE08-64")
print(f"Total citations: {paper.citation_count}")
print(f"Influential: {paper.influential_citation_count}")

cites = get_citations("DOI:10.13182/NSE08-64", limit=50)
for p in cites.papers:
    print(p.summary())  # shows influential count
```

### Build a Bibliography

```python
from tools.research.osti import get_bibtex
from tools.research.crossref import get_work

# BibTeX from OSTI
print(get_bibtex("1893820"))

# Build citation from CrossRef
w = get_work("10.1016/j.anucene.2023.01.001")
print(f"{w.first_author} et al., {w.container_title} {w.volume} ({w.published_date[:4]})")
```

## Rate Limits

| Database | Throttle | Quota |
|----------|----------|-------|
| OSTI | 0.5s | Unlimited |
| arXiv | 3s | Unlimited |
| Scopus | 0.15s | 20,000/week (search) |
| INIS | 0.5s | Unlimited |
| OpenAlex | 0.1s | Unlimited (polite pool with email) |
| CrossRef | 0.1s | Unlimited (polite pool with email) |
| Semantic Scholar | 0.6s | 100 req/min (1000 with key) |
| IAEA-NDS | 0.5s | Unlimited |

All clients enforce throttling automatically.

## Constraints

- **NEVER reference export-controlled codes** (names, manuals, equation numbers)
- **Prefer open-access** when multiple sources cover the same content
- **Verify equations against multiple sources** when possible
- **Map notation** to ORPHEUS conventions when reporting equations
