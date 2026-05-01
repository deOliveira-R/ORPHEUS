---
name: research
description: >
  Use when the user or an agent needs to search scientific literature,
  resolve DOIs, explore citation graphs, or query nuclear data.
  12 databases: OSTI, arXiv, Scopus, INIS, OpenAlex, CrossRef, Semantic
  Scholar, HAL, Zenodo, J-STAGE, IAEA-NDS, EXFOR. Examples: "Find
  papers on collision probability", "Look up Bailey 2009", "What
  cites this DOI?", "U-235 ground state data", "Find Sanchez/Hébert
  HAL deposits", "OpenMC release on Zenodo", "JNST resonance papers",
  "U-235 (n,f) experimental cross sections"
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
| **HAL** | French open archive: CEA/IRSN/CNRS reports, theses, HDRs, OA preprints | None | `tools.research.hal` |
| **Zenodo** | CERN-hosted: code releases (OpenMC/MC-DC/etc.), datasets, OA proceedings — every record has a DOI | None (token optional) | `tools.research.zenodo` |
| **J-STAGE** | Japanese journals incl. JNST (reactor physics, JENDL evals, AESJ proceedings) | None | `tools.research.jstage` |

### Nuclear Data

| Database | Best for | Auth | Module |
|----------|----------|------|--------|
| **IAEA-NDS LiveChart** | Nuclear structure: half-lives, binding energies, gammas, levels | None | `tools.research.iaea_nds` |
| **EXFOR (via IAEA-NDS GitHub)** | Experimental reaction data: measured cross sections (n,f / n,γ / n,2n / ...), angular distributions, fission yields | None (token optional) | `tools.research.exfor` |

## When to Use Which

- **DOE lab report or NUREG** → OSTI
- **Recent/unpublished work** → arXiv
- **Specific journal article** (by DOI, author, journal) → Scopus or CrossRef
- **IAEA reports, Eastern European journals** → INIS
- **CEA, IRSN, CNRS, EDF R&D output; French-authored work; HDR theses** → HAL
- **Open-source nuclear code release** (OpenMC, MC/DC, etc. — version DOI) → Zenodo
- **Cross-section library snapshot, processed JEFF/ENDF dataset, NEA benchmark** → Zenodo
- **PHYSOR / M&C / ICAPP / NURETH proceedings** with author-deposited PDF → Zenodo
- **Japanese reactor-physics work** (Yamamoto, Chiba, Nakagawa) → J-STAGE (JNST)
- **JENDL evaluation papers, AESJ proceedings** → J-STAGE
- **Citation graph** ("what cites X?" / "what does X cite?") → OpenAlex or Semantic Scholar
- **Open-access PDF link** → OpenAlex (`oa_url`), HAL (`pdf_url`), or Zenodo (`files[].download_url`)
- **Impact / influential citations** → Semantic Scholar (`influential_citation_count`)
- **DOI metadata resolution** → CrossRef (`get_work`)
- **Journal-specific search** → CrossRef (`journal_search` by ISSN) or J-STAGE (`issn=`)
- **Broad topic search** → search OSTI + OpenAlex + INIS + HAL + Zenodo in parallel
- **Nuclear structure data** (half-life, binding energy, gammas) → IAEA-NDS
- **Experimental reaction data** (measured σ, dσ/dΩ, fission yields) → EXFOR
- **BibTeX** → OSTI (`get_bibtex`), arXiv (`bibtex_key`), HAL (`get_bibtex`), Zenodo (`get_bibtex`)

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

### HAL (archives-ouvertes.fr)

```python
from tools.research.hal import search, search_simple, get_record, get_bibtex

# Convenience kwargs (recommended starting point)
results = search_simple(author="Alain Hebert", doc_type="ART", year_from=2010)
results = search_simple(author="Hebert")          # single word -> surname
results = search_simple(keywords="MOC pin power", journal="Annals of Nuclear Energy")

# Native Solr query (full power — see field table in hal.py docstring)
results = search('authFullName_t:"Richard Sanchez" AND docType_s:THESE')
results = search('keyword_s:DRAGON5 AND domainAllCode_s:phys.nucl-th', rows=20)

# Restrict to a HAL portal / institution collection
results = search_simple(keywords="criticality", portal="irsn")  # or "cea", "in2p3"

# Resolve by HAL id, DOI, or internal docid
rec = get_record("hal-04810205")
rec = get_record(doi="10.1016/j.anucene.2024.110961")

# Fields populated on every Record
rec.title, rec.authors, rec.doi, rec.journal, rec.volume, rec.pages
rec.publication_year, rec.doc_type, rec.doc_type_label
rec.pdf_url, rec.is_open_access, rec.url, rec.keywords, rec.domains
rec.conference_title, rec.conference_city          # for COMM records
rec.abstract                                       # may be empty

# BibTeX export (uses search API wt=bibtex, NOT the landing page)
print(get_bibtex("hal-04810205"))
```

**docType_s codes** (filter via `doc_type=`): `ART` journal article ·
`COMM` conference paper · `POSTER` · `OUV` book · `COUV` book chapter ·
`REPORT` technical report · `THESE` PhD thesis · `HDR` habilitation ·
`MEM` master's thesis · `PREPRINT` · `SOFTWARE`. Full map in
`tools.research.hal.DOC_TYPES`.

**Field naming**: `_s` = string (exact match, case/diacritic-sensitive),
`_t` = analyzed text (tokenized, diacritic-folded — use for any
user-typed query). Rule of thumb: `_s` for filter fields (HAL ids,
DOIs, controlled codes); `_t` for human strings (titles, names,
journals). The `search_simple` helper picks the right variant
automatically.

**Why HAL** (when other DBs miss it):

- CEA / IRSN technical reports often deposited as `REPORT` here but
  absent from OSTI.
- Many French theses (`THESE`) and HDR monographs (`HDR`) — Sanchez,
  Coste-Delclaux, Reuss, Lautard, etc. — exist nowhere else with
  open-access PDFs.
- Accepted-manuscript PDFs of paywalled Annals of Nuclear Energy
  papers when authors deposit them.

### Zenodo

```python
from tools.research.zenodo import search, get_record, get_bibtex, get_versions

# Free-text + resource-type filter
results = search("OpenMC neutron transport", resource_type="software")
results = search("cross section uncertainty", resource_type="dataset")

# Restrict to a Zenodo community
results = search("PHYSOR 2024", communities="physor")
results = search("benchmark", communities="oecd-nea-data-bank")

# Elasticsearch field syntax for precision
results = search('creators.name:"Romano, Paul" AND title:OpenMC')
results = search('keywords:"covariance" AND publication_date:[2020 TO 2024]')

# Resolve by DOI or numeric record id
rec = get_record(7233164)
rec = get_record(doi="10.5281/zenodo.7233164")

# Per-record fields
rec.title, rec.creators, rec.doi, rec.concept_doi
rec.resource_type, rec.resource_subtype, rec.version
rec.is_open_access, rec.html_url, rec.keywords, rec.communities
rec.files                         # FileEntry(key, size, checksum, download_url)
rec.related_identifiers           # links to papers, source code, datasets
rec.journal_title, rec.conference_title

# Software releases — get every version of a concept DOI
versions = get_versions(7233164)
for v in versions.records:
    print(v.version, v.publication_date)

# BibTeX export (Accept: application/x-bibtex content negotiation)
print(get_bibtex(7233164))
```

**Common resource_type values**: ``"publication"`` · ``"poster"`` ·
``"presentation"`` · ``"dataset"`` · ``"software"`` · ``"image"`` ·
``"video"`` · ``"lesson"`` · ``"physicalobject"`` · ``"workflow"`` ·
``"other"``.

**Auth**: anonymous = 60 req/min (2000/h). Set ``ZENODO_TOKEN`` env
var → 100 req/min (5000/h).

### J-STAGE

```python
from tools.research.jstage import search, JNST_ISSN, JStageAPIError

# JNST = Journal of Nuclear Science and Technology (AESJ).
# Validated 2-parameter combos (the API rejects most 3+ field mixes
# with ERR_001 — see jstage.py docstring for the matrix).
results = search(text="MOC", issn=JNST_ISSN)
results = search(author="Yamamoto", issn=JNST_ISSN, sortflg=2)
results = search(text="resonance self-shielding", pubyearfrom=2018)

# Each Article carries the bilingual JST schema parsed to flat fields
for art in results.articles:
    art.title, art.authors, art.journal_title, art.issn
    art.volume, art.issue, art.start_page, art.publication_year
    art.doi, art.url, art.keywords, art.abstract

# Convenience constants for filter ISSNs
# JNST_ISSN          = "0022-3131"  (print)
# JNST_ISSN_ONLINE   = "1881-1248"
# MEJ_ISSN           = "2187-9745"  Mechanical Engineering Journal
# JPES_ISSN          = "1881-3062"  Journal of Power and Energy Systems

# Catch parameter-combination errors explicitly
try:
    results = search(keyword="MOC", issn=JNST_ISSN)  # ERR_001
except JStageAPIError as e:
    ...
```

**Gotchas**:
- ISSN must include the hyphen: ``"0022-3131"`` not ``"00223131"`` (else ``ERR_008``).
- ``cdjournal`` cannot be used alone — pair with another field.
- ``issn`` cannot combine with ``pubyearfrom``/``pubyearto``.
- Author names should be in romaji (``"Yamamoto"``, ``"Chiba"``); the
  API matches both ``"First Last"`` and ``"Last First"`` orderings.

### EXFOR (Experimental Reaction Data — IAEA-NDS via GitHub)

```python
from tools.research.exfor import (
    list_datasets, get_dataset, get_entry_json, parse_filename,
    COMMON_REACTIONS, MT_MAP,
)

# 1. Discover every measurement of a reaction
files = list_datasets(projectile="n", target="U-235",
                      reaction="(n,f)", observable="xs")
# Returns list[DatasetRef] — one per EXFOR subentry.
# Loose target/reaction forms are normalised:
#   "U235" / "u-235" / "235U" -> "U-235"
#   "(n,f)" / "n,f" / "N,F" -> "n-f"

# 2. Fetch and parse the tabulated data
ds = get_dataset(files[0])
ds.x, ds.dx, ds.y, ds.dy           # parallel float lists
ds.column_labels                   # e.g. ['E_in(MeV)', 'dE_in(MeV)', 'XS(B)', 'dXS(B)']
ds.target, ds.process              # '92-U-235', 'N,F'
ds.mf, ds.mt                       # ENDF MF/MT (e.g. '3', '18')
ds.author, ds.year, ds.facility, ds.institute
ds.entry_id, ds.reaction_code, ds.reference
ds.master_file_url, ds.nds_url     # canonical EXFOR sources

# 3. Pull the full ENTRY JSON (BIB block + every subentry)
entry = get_entry_json("14015")
entry["bib_record"]["AUTHOR"]
entry["data_tables"]               # all subentries with original units

# Filename parser (also called automatically by list_datasets)
parse_filename("U-235_n-f_A.D.Carlson-14015-002-0-1991.txt")
# -> {'target': 'U-235', 'reaction': 'n-f', 'author': 'A.D.Carlson',
#     'entry_id': '14015', 'subentry': '002', 'pointer': '0', 'year': '1991'}
```

**Reaction directory codes** (``COMMON_REACTIONS`` maps friendly →
directory): ``n-f`` fission · ``n-2n``, ``n-3n``, ``n-4n`` · ``n-g``
capture · ``n-tot`` total · ``n-el`` elastic · ``n-inl`` inelastic ·
``n-non`` nonelastic · ``n-abs`` absorption · ``n-p``, ``n-a`` ·
``n-x`` total emission · ``n-sct`` total scattering.

**Observables**: ``xs`` (cross section, default) · ``angle``
(angular distribution) · ``energy`` (secondary spectrum) ·
``fission/yield/cumulative`` · ``fission/yield/independent`` ·
``fission/yield/primary`` · ``neutron`` (selected neutron observables).

**Data sources** (both maintained by IAEA-NDS, MIT licensed):
- ``IAEA-NDS/exfortables_py`` — pre-tabulated ``(x, dx, y, dy)`` files.
- ``IAEA-NDS/exfor_json`` — full ENTRY/SUBENTRY JSON dumps.

GitHub raw content is unrestricted; the GitHub *contents* API
(used by ``list_datasets``) is throttled at 60/h anonymous, 5000/h
with ``GITHUB_TOKEN`` set.

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
   - HAL: French OA + theses    Zenodo: code/data/proceedings DOIs
   - J-STAGE: JNST + AESJ       (CEA/IRSN/CNRS reports & HDRs in HAL)
3. Deduplicate by DOI (and HAL halId / Zenodo record id where DOI is missing)
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
| HAL | 0.5s | Unlimited (Solr API; default 30 rows, max 10000) |
| Zenodo | 1.1s | 60/min anonymous (2000/h); 100/min with `ZENODO_TOKEN` |
| J-STAGE | 0.5s | Unlimited (max 1000 results per page, ATOM XML only) |
| IAEA-NDS | 0.5s | Unlimited |
| EXFOR (GitHub) | 1.0s | Raw downloads unlimited; contents API 60/h anon, 5000/h with `GITHUB_TOKEN` |

All clients enforce throttling automatically.

## Constraints

- **NEVER reference export-controlled codes** (names, manuals, equation numbers)
- **Prefer open-access** when multiple sources cover the same content
- **Verify equations against multiple sources** when possible
- **Map notation** to ORPHEUS conventions when reporting equations
