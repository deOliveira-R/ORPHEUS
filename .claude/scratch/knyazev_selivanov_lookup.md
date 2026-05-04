# Knyazev-Selivanov 2014 — Bibliographic Lookup

## Headline finding

**The "Knyazev-Selivanov 2014" reference does not exist as cited.** What
exists — and what every load-bearing claim in the ORPHEUS memory chain
actually rests on — is a single-authored paper by **A. P. Knyazev** in
**Atomic Energy, 1993**. The "2014" date and "Selivanov" co-author appear
to be ghosts introduced by a downstream sub-agent's wish-list memo and
propagated forward; they have no basis in any database (CrossRef,
OpenAlex, Semantic Scholar, OSTI, INIS, web search) or in any of the
in-tree memory files that originally cited the work.

**The canonical reference the user actually needs:**

> **A. P. Knyazev (1993).** "Solution of the transport equation in
> integral form in a one-dimensional cylindrical geometry with linearly
> anisotropic scattering." *Atomic Energy* (Plenum/Springer translation
> of *Atomnaya Énergiya*), **Vol. 74, No. 5, pp. 368–374**, May 1993.
> DOI: **`10.1007/BF00844623`**.
>
> Russian original: *Atomnaya Énergiya*, Vol. 74, No. 5, pp. **403-410**,
> May 1993. Translation note from Springer's Additional Information
> section: "Institute of Technical and Experimental Physics. Translated
> from *Atomnaya Énergiya*, Vol. 74, No. 5, pp. 403-410, May, 1993."

## Verified citation

| Field             | Value                                                                                                                         |
| ----------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| Authors           | **A. P. Knyazev** (sole author; not Knyazev *and* Selivanov)                                                                  |
| Affiliation       | Institute of Technical and Experimental Physics                                                                               |
| Title (English)   | Solution of the transport equation in integral form in a one-dimensional cylindrical geometry with linearly anisotropic scattering |
| Journal (English) | Atomic Energy (Plenum / Springer)                                                                                             |
| Volume / issue    | 74 / 5                                                                                                                        |
| Pages (English)   | 368-374                                                                                                                       |
| Year / month      | 1993, May                                                                                                                     |
| DOI               | `10.1007/BF00844623`                                                                                                          |
| Russian original  | *Atomnaya Énergiya*, Vol. 74, No. 5, pp. 403-410, May 1993                                                                    |
| ISSN (English)    | 1063-4258 (print)                                                                                                             |

## DOI verification

CrossRef metadata API (`https://api.crossref.org/works/10.1007/BF00844623`)
returned the full record: title, single-author (A. P. Knyazev), journal
*Atomic Energy*, vol. 74 issue 5, pp. 368-374, published 1993-05.

Springer landing page resolves correctly:
**https://link.springer.com/article/10.1007/BF00844623**

The `Cite this article` block on the Springer page reads:

> Knyazev, A.P. Solution of the transport equation in integral form in a
> one-dimensional cylindrical geometry with linearly anisotropic
> scattering.

confirming sole authorship.

## Open access status

**Not open access.** Springer paywall. No OA mirror found. Semantic Scholar
reports `is_open_access: False`, `oa_pdf_url: None`. OpenAlex returns
`is_open_access: False`, `oa_url: ''`. No arXiv preprint. No HAL or Zenodo
deposit.

The user MAY have institutional access to Springer/Atomic Energy via a
nuclear-engineering library subscription. The Russian original
(*Atomnaya Énergiya* 74, 403-410) MAY be locatable via a Russian
academic library or a Soviet-era reactor-physics archive, but that's
a Russian-language source.

## What the paper provides

Based on the title and the surrounding citation context in
`/workspaces/ORPHEUS/.claude/agent-memory/numerics-investigator/issue_132_cylinder_hebert.md`,
the paper provides the **closed-form integral expression for the
linearly-anisotropic-scattering 1-D cylinder collision/transmission
kernel**, expressed via Bickley-Naylor higher-order functions
**Ki_{2+k}**. This is the polar-angle-integrated chord kernel that ORPHEUS
needs as a 3-D correction to the cylinder Hébert closure: the cylinder
P-block currently uses Ki_3 (the Sanchez 1982 §IV.A canonical form),
but the BC kernel `compute_G_bc` uses a 2-D-projected-cosine factor
that under-predicts row-sum partition by ~7.6%. Knyazev 1993 provides
the analytic Ki_{2+k} machinery that closes this gap by integrating
out θ_p analytically against the µ-weighted partial-current measure.

Per the existing numerics-investigator memo:

> Polar integration of `∫_0^π sin θ_p · µ_3D · P̃_m(µ_3D) P̃_n(µ_3D) ·
> e^{-τ_2D / sin θ_p} dθ_p` against the µ-weighted partial-current
> measure gives the Knyazev Ki_(3+k_m+k_n) expansion (one Ki order
> HIGHER than the corresponding cylinder P/G primitives, due to the
> extra µ_3D = sin θ_p factor).

So the deliverable inside the Knyazev 1993 paper is, structurally,
the cylinder analog of Sanchez 1986's sphere Eq. (A6) — a 1-D
quadrature evaluator for the homogeneous linearly-anisotropic
cylinder kernel, written in Bickley-Naylor functions of order 2+k.
This is exactly the missing primitive cited as a wish in the gap
analysis.

## Where I found it

- **CrossRef** (`tools.research.crossref.get_work("10.1007/BF00844623")`)
  returned the full canonical record. **Authoritative source for this
  citation.**
- **OpenAlex** (`tools.research.openalex.get_work("10.1007/BF00844623")`)
  confirmed authorship and OA status.
- **Springer landing page** (HTML scrape of
  `https://link.springer.com/article/10.1007/BF00844623`) provided the
  Russian-original journal/volume/page mapping plus the institutional
  affiliation — neither of which is in the CrossRef metadata.
- **In-tree references** (already known to the user but worth flagging
  for traceability):
  - `/workspaces/ORPHEUS/.claude/agent-memory/numerics-investigator/issue_132_cylinder_hebert.md` —
    cites `Knyazev 1993, Atomic Energy 74, DOI 10.1007/BF00844623`
    correctly. **This is the original, accurate citation.**
  - `/workspaces/ORPHEUS/.claude/agent-memory/numerics-investigator/_archive/specular_mb_phase4_cyl_slab.md` —
    cites the same Knyazev Ki_{3+k} machinery without a year, but
    consistent with the 1993 paper.
  - `/workspaces/ORPHEUS/.claude/agent-memory/literature-researcher/peierls_greens_extensions_lit.md` —
    introduces the **incorrect** "Knyazev, B. A. & Selivanov, A. N. (2014)"
    Bickley-Naylor functions wish-list item. This entry is the source
    of the phantom 2014 citation.
  - `/workspaces/ORPHEUS/.claude/agent-memory/literature-researcher/sanchez_chandrasekhar_three_meanings.md` and
    `/workspaces/ORPHEUS/.claude/scratch/sanchez_chandrasekhar_gap.md` —
    propagated the phantom 2014 citation forward, sometimes alongside
    the legitimate 1993 paper, creating the appearance of two distinct
    works.

## Confidence

**High** for the 1993 paper:
- DOI resolves on Springer.
- CrossRef and OpenAlex both return the full record under the same DOI
  with single authorship by A. P. Knyazev.
- The Springer landing page exposes the Russian original
  (*Atomnaya Énergiya* 74, 403-410, May 1993) and the institutional
  affiliation.
- The 1993 citation appears in `issue_132_cylinder_hebert.md` with a
  load-bearing technical claim (the Ki_{2+k} structure is what the
  cylinder Hébert closure requires) and is internally consistent with
  Sanchez 1982 §IV.A which the same memo cross-references for the
  Ki_1 + Ki_3 structure.

**High** for the negative finding (no Knyazev-Selivanov 2014 paper):
- CrossRef searches over multiple author/keyword combinations return
  zero matches.
- OpenAlex returns zero hits even on the broader
  `Knyazev Selivanov` author pair.
- Semantic Scholar returns zero relevant hits.
- INIS rejected the connection (transient SSL errors mid-search, but
  the alternate routes all returned zero).
- OSTI returned irrelevant results (Selivanov is a common Russian
  surname; matches were either physicists in unrelated fields or
  accelerator-conference proceedings).
- Web search ("Knyazev Selivanov" + Bickley + Bickley-Naylor + Nauka +
  Moscow + variations on the year) returned zero hits across regular
  Web, Google Scholar shadow indexes, and Russian-academic shadow
  indexes.
- The journal-restricted query against `Atomic Energy` (ISSN 1063-4258)
  via CrossRef returned **only one** A. P. Knyazev paper — the 1993
  one. There is no 2014 follow-up, no Selivanov co-authorship, and no
  Knyazev-Selivanov collaboration anywhere in the *Atomic Energy*
  publication record.

The "2014" wish-list entry in `peierls_greens_extensions_lit.md`
describes a *book* — "Bickley-Naylor functions and their integrals,
Nauka, Moscow / English translation" — that does not appear in any
catalogued bibliographic source. It is most likely a hallucinated
extension of the legitimate Knyazev 1993 paper into a fictional
co-authored monograph. Russian-academic publisher Nauka does exist;
no such 2014 book by either author is catalogued in WorldCat,
COPAC/Library Hub, or the Russian State Library OPAC indices that the
above search infrastructure routes through.

## Recommendation

The user should:

1. **Use the 1993 paper as the canonical citation.** It is the actual
   technical source. The DOI `10.1007/BF00844623` resolves cleanly.
2. **Remove or annotate the "Knyazev-Selivanov 2014" entries** in
   `/workspaces/ORPHEUS/.claude/agent-memory/literature-researcher/peierls_greens_extensions_lit.md`,
   `/workspaces/ORPHEUS/.claude/agent-memory/literature-researcher/sanchez_chandrasekhar_three_meanings.md`,
   and `/workspaces/ORPHEUS/.claude/scratch/sanchez_chandrasekhar_gap.md`.
   Each of those files should be updated to point at the 1993 DOI.
3. **Acquire the 1993 paper from Springer / institutional access**, or
   request the Russian original (*Atomnaya Énergiya* 74, 403-410)
   via a library that holds Soviet-era nuclear-physics journals.
   Both are paywalled; the English Springer version is the more
   accessible.

If the user has independent reason to believe a different
Knyazev-Selivanov publication exists (e.g., a citation snippet from
another paper they have in hand), I should be re-tasked with that
specific source citation — bibliographic information from a citing
paper would over-rule the negative finding above.

## Alternatives if the user wants a different cylinder Meaning-2 source

Beyond Knyazev 1993 (which is a closed-form 1-D linearly-anisotropic
cylinder integrator) and Carlvik 1965 (which the user reports they
already have — the Wigner-Seitz CP foundational paper), the candidate
cylinder Meaning-2 references are:

1. **Sanchez 1982** *NSE* 80 (cited in `issue_132_cylinder_hebert.md`
   §IV.A for the Ki_1 + Ki_3 canonical structure) — a textbook-level
   walk-through of cylinder CP via Bickley-Naylor primitives. Likely
   covers the same Ki_{2+k} structure as Knyazev 1993 but in less
   focused form. Verify with CrossRef before relying.
2. **Hébert 2009** book chapters 6-9 — the modern reactor-physics
   monograph treatment. Not local in full; only Ch.3 is in
   `scratch/literature/`. This is in the existing wish-list.
3. **Garcia 2006** *Annals of Nuclear Energy* — "Collision
   Probabilities in r-θ-z Geometry." This is **LOCAL** in
   `scratch/literature/Garcia(2006)Collision Probabilities in r-theta-z Geometry.pdf`.
   3-D cylindrical geometry CP, more general than 1-D infinite
   cylinder. Worth re-reading if the goal is cylinder Meaning-2.
4. **Knyazev 1993 itself** is the most direct match — the only one
   with the explicit linearly-anisotropic + Ki_{2+k} structure.

The cylinder Meaning-2 acquisition gap is therefore reduced from
**two unknown references** (Carlvik 1965 + Knyazev-Selivanov 2014) to
**one known reference** (Knyazev 1993, DOI `10.1007/BF00844623`,
paywalled but accessible via Springer institutional subscription).
