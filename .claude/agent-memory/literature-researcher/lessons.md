# Literature Researcher — Lessons

Behavioural corrections only: a mistake, and the check that now prevents it.
The definition (`docs/development/agents/literature-researcher.md`) and the
`research` skill hold the source priorities, the notation mapping and the
extraction procedure; nothing here repeats them. Records of what a paper says
live in the topic files that `MEMORY.md` indexes.

Sixteen lessons, L-001 to L-016, once stood here. An
audit on 2026-09-22 retired all sixteen as carried by the definition, the skill
or a rule. A re-check on 2026-11-16 read each lesson against those carriers
(the pre-retirement text is `git show 21f999c0:.claude/agent-memory/literature-researcher/lessons.md`):
eight are carried whole and stay retired, and eight had a part that no carrier
states. Those eight parts are below under their original numbers. The ledger at
the end names the carrier of everything else. The war stories are in git
history and in the topic files named. A new lesson takes L-017 or later.

Markers: `[M]` measured, with its date; `[R]` reasoned or copied forward, not
re-run on 2026-11-16 (that pass had no web tool).

---

## L-002 — Enumerate the journal slot; do not only search the title

The definition says to check a citation by its slot (journal, volume, first
page). This is how the slot is read.

- check: (a) list the venue's volume through CrossRef, filtered by ISSN and
  publication date and selecting title, volume, page and DOI, then filter the
  result on the first page; the record at that page is the slot's true
  occupant `[R]`. (b) Get the author's count from the journal index,
  `journals/<ISSN>/works?query.author=<surname>`, whose `total-results` field
  is a denominator that can be reported and that a homograph cannot inflate
  `[R]`.
- worked case: `bailey_ane_chimera_citation_refuted.md`.

## L-003 — A section or equation pointer inside a correct citation is its own claim

A brief gave the right book and chapter with the wrong subsection: it cited
the spherical section's equations for the cylinder. The book was right, so the
citation looked verified.

- check: when a brief hands over a section and an equation number together,
  search the OCR sidecar, in Python, for the section headings and for the
  equation number before opening any page, and confirm that the equation falls
  inside the section named and that the section treats the geometry named.
  A sidecar writes an equation number either as `\tag{<number>}` or as a
  trailing `(<number>)` inside the display (`[M]` 2026-11-16: 72 of 77
  sidecars use `\tag{}`; the Reed-Lathrop 1970 sidecar uses the other form
  and has 0 `\tag{` matches). Search both forms, and first search for an
  equation number known to be in the paper: if that finds nothing, the pattern
  is wrong, not the pointer.

## L-009 — Read the open-access status value, and follow a moved archive

The definition says that a publisher's `403` is usually a bot challenge and
that the open-access status decides. The status values mean different user
actions.

- check `[R]` (field names copied forward, not queried again on 2026-11-16):
  `bronze` means the publisher's own site serves the paper free; tell
  the user to download it in a browser. Its sign in OpenAlex is an `oa_url`
  that carries `?needAccess=true`, which is not a false flag. `green` means a
  repository copy exists and `locations[]` names it. A record with
  `any_repository_has_fulltext: false` and a single location has no repository
  copy, so stop searching repositories.
- check: a `404` at the journal's expected host can mean that the archive
  moved to another publisher. Follow the DOI before reporting the paper absent.
- worked case: `ty_polar_quadrature_moc.md`.

## L-010 — A fourth check on a printed equation: is it vacuous?

The definition lists three checks beyond the scan (the twin equation, the seed
or terminus, a reproduced table). A fourth costs nothing.

- check: read the equation for vacuity. A recursion whose two sides name the
  same symbol, a definition that defines nothing and a sum that cannot close
  are printing errors, whatever the scan shows.

## L-011 — When a symbol's meaning is ambiguous, the paper's table decides

Two readings of a printed sum were typographically identical, and one of them
made the quantity vanish for every scheme. Reading the page again could not
choose between them.

- check: when an extraction depends on which quantity a symbol denotes, or
  when two sources seem to contradict each other, code each candidate reading
  and reproduce a table the paper prints. The reading that reproduces the
  printed digits is the paper's. Report the reproduction as `[M]`.
- check: before recommending a paper's diagnostic as a test, compute it on the
  configuration the user actually runs. In the worked case the diagnostic was
  round-off for every convention on the user's quadrature, so it would have
  been a test that cannot fail.
- worked case: `morel_montry_tau_angular_cell_edges.md`.

## L-012 — A "no source does X" finding: who is credited, what mechanism, and how old is the search

The count of sources behind a negative is the instrument-doctrine rule's
statement X2, and the definition covers a search tool that cannot see the
field. Three checks are carried by neither.

- check: search the local text for the author it credits with the
  construction ("the approach proposed by …"). That credited primary is
  usually not among the sources read. Obtain it, or say in the finding that it
  is unread.
- check: restate the negative as a mechanism. "All three sources define the
  face cosine by a conservation recursion, and whether that lands on the
  geometric value depends on the quadrature" can be checked and reused; "no
  source uses the geometric value" cannot.
- check: a negative is as old as the library it was searched in. When
  `scratch/literature/` has gained files since the negative was written,
  search again, and include every sidecar suffix: a sidecar made from a PDF's
  own text layer is named `<stem>.textlayer.md`, and a search limited to the
  OCR stems misses it (`[M]` 2026-11-16: 2 such sidecars among 77).
- worked case: `issue_158_linear_discontinuous_cell_update.md`.

## L-013 — Where authors state what their scheme costs, and how to test "X determines Y"

The definition says to count a named scheme's equations and unknowns. Two
further checks belong to the same reading.

- check: the authors' own statement of what their choice costs, and of which
  geometries they did not treat, is in the Conclusions (often as future work)
  or in the Introduction's paragraph that positions the work against its
  rival. Read both places before reporting what a paper covers.
- check: a later paper often compresses an earlier derivation into "condition
  X determines Y". Count the equations in X against the unknowns in Y. If X
  has fewer, exhibit a second solution numerically, report the claim as
  necessary and not sufficient, and find the earlier paper's actual
  determining condition.
- worked cases: `reed_lathrop_1970_angular_truncation.md`,
  `morel_montry_tau_angular_cell_edges.md`.

## L-016 — "Is my measurement a property of the format, the evaluation or the processing?" is three questions

A processing code's manual stopped short of the answer: it named the
subroutine that chooses a default and did not print the subroutine's table.

- check: answer each of the three from its own source: the format's manual
  (what is permitted), the evaluated file (what this nuclide contains) and the
  processing code's manual together with the processed file (what was done to
  it). Report where the three disagree.
- check: when the manual under-determines the behaviour, parse the real
  evaluated files behind the user's data. A file in the user's own set that
  behaves differently under the same processing run is the strongest evidence
  available.
- check: reproduce the user's headline number from the raw upstream data, and
  state for how many of the user's cases the confirmation holds.
- worked case: `endf6_gendf_njoy_n2n_formats.md`.

---

## Ledger — where the retired lessons are carried (re-checked 2026-11-16)

"Definition" is `.claude/agents/literature-researcher/AGENT.md`, generated from
the docs page named above; "skill" is `.claude/skills/research/SKILL.md`. This
agent loads no rule, so a rule carries a lesson only in a dispatch whose brief
lists that rule.

| Lesson | Subject | Carrier |
|---|---|---|
| L-001 | name the geometry when mapping a direction cosine; the cylinder map | definition §3 |
| L-002 | a citation known only from a memo is unverified; the fused citation | definition §2 "Provenance"; skill "Phantom-citation guard"; the slot recipe is kept above |
| L-003 | a docstring's citation is a claim | definition §2 "Provenance"; the section pointer is kept above |
| L-004 | a results catalogue is not a method source | definition §2 "Method, not topic"; skill "Catalogue ≠ method-source" |
| L-005 | classify a method, and a geometry, from the paper's own equations | definition §2 "Method, not topic"; the future-work clause is kept in L-013 |
| L-006 | the local library first; never substitute a source unasked | definition §1; the workflows rule, W7 |
| L-007 | a dead Zotero server | definition §1 item 2 |
| L-008 | table bodies are in the `.mocr.json` | definition §1 item 1 |
| L-009 | a `403` is not a paywall | definition §1 "Not in the local library"; the status values are kept above |
| L-010 | the print itself can be wrong | definition §2 "The equation"; the vacuity check is kept above |
| L-011 | reproduce the paper's table | definition §2 "The equation", as a third check; the trigger and the corollary are kept above |
| L-012 | a negative needs its count; a search tool blind to the field; the citation graph | instrument-doctrine X2 (when briefed); definition §1 "Searching"; definition "Return" for the sites a refuted claim reached; three checks are kept above |
| L-013 | count a named scheme's equations and unknowns | definition §2 "A named scheme"; two checks are kept above |
| L-014 | two sources conflict: read the later one's citation of the earlier | definition §2 "A disagreement between two sources" |
| L-015 | a format question is answered from the standard; search extracted text in Python; layout loss in fractions | definition, opening paragraph, §1 and §2 "The equation" |
| L-016 | measure the data when the manual stops | not carried; kept above |
