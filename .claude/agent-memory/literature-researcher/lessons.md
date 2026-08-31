# Literature Researcher — Lessons

Behavioral corrections only (failure → rule). The AGENT.md / `research`
skill hold source priorities, notation mapping, and extraction
procedure — never duplicate those here. Durable paper-extraction
RECORDS (equations, notation maps, value tables) live in the topic
files indexed by `MEMORY.md`, NOT here.

---

## L-001 — Bailey is the canonical cylindrical SN reference; state the geometry when mapping µ

Bailey, Morel & Chang **(2010)** NSE **165**(2):149-169,
`10.13182/NSE08-66` (LOCAL): §V is R-Z, and **Eq. 50** = α dome
recursion, **Eq. 52** = the per-ξ-level edge-cosine recursion (⟸ the
η-ascending level ordering lives HERE, not in 50), **Eq. 74** =
Morel-Montry τ. Sphere twins = Eqs. 11/12.

**Why:** Bailey's µ = radial cosine = ORPHEUS `mu_x`, but in
slab/spherical contexts µ means the *streaming* cosine. The symbol is
overloaded across geometries.

**How to apply:** ALWAYS name the geometry when mapping µ. For
curvilinear SN balance / the ΔA/w geometry factor, start at Bailey §4,
not Lewis & Miller (who omit the asymptotic analysis explaining *why*
the geometry factor exists).

⛔ **CORRECTION (2026-08-27) — AGENT.md's own notation table is INVERTED
for the cylinder; L-001 above is the correct one.** `[VERIFIED ON SCAN]`
Hébert 2009 Ch.3 (3.152)/(3.157) p.91 and BMC 2010 Eq. (48) p.156 use the
**SAME** assignment: **µ = RADIAL, η = AZIMUTHAL, ξ = AXIAL (level)**.
AGENT.md says "our `mu_x` = their η (radial); our `mu_z` = their µ
(axial)" — backwards for BOTH named sources. The redistribution term is
the tell and it is one grep: whichever symbol sits inside `∂/∂ω(·ψ)` is
the **azimuthal** cosine (η in both), and whichever multiplies
`∂/∂r(rψ)` is the **radial** one (µ in both). Run that grep before
mapping any curvilinear symbol; do not trust the inherited table.
See [[bailey-ane-chimera-citation-refuted]].

## L-002 — A citation that lives only in a memory is a PHANTOM; resolve every DOI to a real database

→ standing directive now in AGENT.md §6 ("Verify provenance"). War-story kept for forensic value.

"Knyazev-Selivanov 2014" was hallucinated in a wish-list memo, then
re-cited forward across later memos as if real. It exists in NO
database. The real source is sole-author **A.P. Knyazev 1993** *Atomic
Energy* 74(5):368-374, DOI 10.1007/BF00844623.

**Why:** a self-cite inside the memory graph manufactures false
authority — each forward reference makes the phantom look more
established. The error compounds silently.

**How to apply:** before treating ANY cited paper as real, resolve it
against CrossRef/OpenAlex/OSTI. A reference whose only provenance is a
prior memory entry is UNVERIFIED until a database confirms author +
year + DOI. See [[knyazev-1993-cylinder-anisotropic-ic]].

⭐ **Sharpening (2026-08-27) — the CHIMERA, which defeats a title search
because every field is REAL, just from a different paper. Test the
VOLUME/PAGE SLOT, not the title.** `[M]` ORPHEUS cited "Bailey, Adams,
Yang, Zika (2009), *A piecewise linear discontinuous FE spatial
discretization of the transport equation*, Ann. Nucl. Energy **35**,
1929-1936". Title search is inconclusive-looking (a near-title paper DOES
exist — an LLNL conference paper, different authors, geometry clause
dropped); author search finds a real quartet (from a *different*,
already-retracted JCP diffusion paper); the page range makes it look
*more* precise than an ordinary phantom. **The decisive test is the
slot**: a journal+volume+pages triple is a coordinate that is either
occupied by your paper or by someone else's. `[M]` ANE 35:1929-1936 is
Zio & Zoia's Bayesian-MCMC BWR paper (`10.1016/j.anucene.2008.03.007`).

**Three mechanical steps, all cheap:**
(a) **Enumerate the venue's slot** — CrossRef
`filter=issn:<ISSN>,from-pub-date:…` with `select=title,volume,page,DOI`,
then filter locally on the first page. Names the true occupant.
(b) **Check volume↔year internally** — a volume number and a year are two
claims about the same object; ANE 35 = **2008**, so "35 … (2009)" was
self-refuting before any lookup.
(c) **Get the denominator from the journal index, not free text** —
`journals/<ISSN>/works?query.author=<surname>` returns `total-results`
(here **4**, across the entire run of ANE, none a match). That is a
countable, reportable denominator (L-012a) and it is immune to the
homograph blindness of free-text search (L-012 sharpening).

⭐ And the payoff that makes this worth doing rather than just deleting the
citation: **a chimera's EQUATION NUMBERS are usually still valid — for the
paper the number was stolen from.** Here "Eq. 50" was correct all along for
BMC 2010 NSE 165, which was LOCAL and already in `refs.bib`. So the fix was
a **re-point, not a deletion**, and reporting it as "citation is fake,
delete it" would have thrown away correct physics provenance. Always ask
*which real paper does this equation number fit?* before recommending
removal. See [[bailey-ane-chimera-citation-refuted]].

## L-003 — A code docstring's citation is a CLAIM, not a fact; verify it

→ standing directive now in AGENT.md §6 ("Verify provenance"). War-story kept for forensic value.

`reduced_operator.py` cited Adams-Yang-Zika 2008 JCP (a *diffusion*
paper) for what is actually Bailey-Morel-Chang 2010 NSE 165 *angular*
SN machinery. Wrong paper, wrong sub-field.

**Why:** docstrings drift from the math they annotate; an
authoritative-looking inline citation can point at the wrong work and
mislead a whole investigation. See [[sphere-sn-pole-closure-canonical]].

**How to apply:** when a code comment names a reference, confirm the
equation actually lives in THAT paper before reusing the citation in
output. Flag the mismatch to the user.

⭐ **Sharpening (2026-08-11): the SECTION/EQUATION pointer inside an
otherwise-correct citation is its own claim, and it fails independently.**
A brief said "Hébert §3.9.4, pp. 141-144, Eqs. 3.418-3.439" for the
**cylindrical** closure — right book, right chapter, wrong subsection:
§3.9.4 is the SPHERE; the cylinder is **§3.9.3** and the whole
3.418-3.439 range is spherical. The same brief called (3.437)/(3.439)
"a weighted one" when they are (3.431) rearranged. Both were caught by
**two greps before any page was opened**: `grep -n "^#.*<sec>"` on the
sidecar for the section-header index, then `grep -n "tag{<eq>}"` to see
which section each equation actually falls in. Do this FIRST on any
brief that hands you a section+equation pair — an inherited pointer sends
you to the wrong geometry's algebra, which is the Cardinal-Rule-2
failure the extraction exists to prevent.

## L-004 — A results CATALOGUE is a test set, not a method source; the derivation is in its cited literature

→ standing directive now in AGENT.md §6 ("Classify by the source body"). War-story kept for forensic value.

Sood-Forster-Parsons (LA-13511 1999 / Prog. Nucl. Energy 2003) is a
75-problem benchmark catalogue + transport-equation definitions. It
contains NO F_N / Wiener-Hopf / X-function / multi-region / anisotropic
derivations — those live in the *cited* primary sources (KLL 1974,
Grandjean-Siewert 1979, Siewert-Thomas 1986, Westfall-Metcalf,
Burkart-Ishiguro-Siewert, Atalay). The 2003 journal port added nothing
methodological beyond a general-multigroup k_inf appendix.

**Why:** the user reasonably hypothesised the journal version "hid" the
method; reading both confirmed it does not. Hunting a derivation inside
a catalogue burns time.

**How to apply:** use a benchmark report for VALUES and the reference
graph; chase METHOD to the primary papers it cites. See
[[sood-2003-vs-1999-extraction]] and [[sood-fn-method-full-extraction]].

## L-005 — Classify a paper's METHOD by reading it, not by the context that cited it

→ standing directive now in AGENT.md §6 ("Classify by the source body"). War-story kept for forensic value.

Atalay 1997, Burkart-Ishiguro-Siewert 1976, and Dahl-Sjöstrand 1979 all
arrive via Sood's "F_N" reference cluster, yet NONE is F_N: Atalay/BIS
are **Case singular-eigenfunction** (+ Fredholm / Chandrasekhar
H-function), Dahl-Sjöstrand is **Legendre-Galerkin on Carlvik's
integral equation**. This is load-bearing for ORPHEUS architecture —
they do NOT extend `fn_method/core/` (different machinery: Case
eigenfunctions / Galerkin vs Marshak-moment collocation).

**Why:** mis-routing a Case-method paper into an F_N package grafts the
wrong mathematical machinery onto the codebase (Cardinal Rule 2). The
citing context predicts topic, NOT method.

**How to apply:** before recommending where a reference lands in the
tree, state the actual solution method from the paper body. F_N =
boundary collocation on Marshak moments; Case = singular
eigenfunctions; Galerkin = polynomial expansion of the integral form —
they are structurally independent and that independence is also a V&V
asset (cross-method cross-checks). See [[atalay-1997-reflected-anisotropic]],
[[burkart-ishiguro-siewert-1976-two-region-anisotropic]],
[[dahl-sjostrand-1979-anisotropic-slab-sphere]].

⭐ **Sharpening (2026-08-25) — the GEOMETRY is a method property too, and a
title + a citing context can agree with each other and both be wrong.** I
reported Morel-Wareing-Smith 1996 (JCP 128:445, *"A Linear-Discontinuous
Spatial Differencing Scheme for S_n Radiative Transfer Calculations"*) as
**the 1-D SPHERICAL spatial-LD primary** without the PDF. Two things
licensed it and neither was evidence: the dispatch brief hypothesised
"1-D spherical", and Lathrop 2000 cites it as ref 1 for *"improving
spatial difference representations"* — a claim about TOPIC, not geometry.
`[M]` on the acquired PDF: its Eq. (53) is `μ ∂ψ/∂z` (**slab**), it has
**no angular-derivative term anywhere**, its abstract says "1D
slab-geometry" three times, and its **§8 names 1-D spherical geometry as
FUTURE WORK** — the exact opposite of the claim. It shipped into a verdict
line and a memory topic file.

**Why this is worse than the L-005 base case:** a wrong METHOD family
usually shows up as soon as you open the paper looking for the machinery.
A wrong GEOMETRY does not — the method name ("linear-discontinuous"), the
author set, the journal and the era all matched perfectly, and the paper
really is the lumping primary the brief wanted. Only the coordinate system
was wrong, and nothing in the citation graph encodes it.

**How to apply:** when a brief or a citing paper supplies the GEOMETRY,
treat it as an unverified hypothesis like any other. The check is one grep
and it is decisive: **grep the paper's own transport equation for the
streaming operator** (`∂/∂z` / `∂/∂x` vs `∂(r²ψ)/∂r` / `∂(rψ)/∂r`) and for
the presence of an **angular-derivative term**. A paper with no
`∂/∂μ`/`∂/∂ω` term cannot be curvilinear, whatever its title says. Then
read the CONCLUSIONS for a "future work" clause — that is where an author
states which geometries they did NOT do (L-013(c): the concession lives in
the conclusions or the intro, never the body).

## L-013 — A NAMED scheme is ambiguous until you count its equations: the same closure can differ in which symbols are UNKNOWN

`[M]` 2026-08-11, Reed & Lathrop 1970 vs Morel-Montry/BMC. Both are called "the
weighted diamond". Both print the **identical** closure
`ψ_m = τ ψ_{m+1/2} + (1−τ)ψ_{m−1/2}` and the **identical** barycentric node
definition, both credit **Grant 1968**, and both use the letter τ for the same
object. Reading either equation confirms nothing. The schemes differ because R&L
impose a **third** equation (the angular consistency condition, Eq. 13b) that
Morel-Montry does not — which flips the **ordinates from inputs to outputs**: R&L
take the weights and SOLVE for `μ_m` (Eq. 14, a quadratic per ordinate), so the
quadrature is a product of the difference scheme. The brief invited exactly the
wrong answer ("what is R&L's weighted diamond?" beside "ORPHEUS implements the
Morel-Montry weighted diamond"), and reporting them as one scheme would have
grafted the wrong provenance onto a shipped closure.

**Why:** L-005 discriminates method FAMILIES (F_N vs Case vs Galerkin) and is
answered by reading what machinery the paper builds. This is the sub-family case,
where the machinery, the symbols, and the printed equations all MATCH and the only
difference is the **unknown set**. No amount of equation-matching detects it; only
counting the constraints and asking "solved for what?" does.

**How to apply:** for any scheme the brief names, (a) **enumerate every equation
the paper imposes** and (b) **state which symbols are unknowns and which are
given** — write both into the deliverable as a table against the local
implementation. Then (c) **read the authors' own statement of what their choice
costs** — a scheme that constrains extra unknowns always pays somewhere, and the
primary usually admits it plainly where a later critic states it as a dig. It
lives in one of exactly two places: the **CONCLUSIONS**, or the **INTRODUCTION's
positioning-against-the-rival paragraph** — never the body. (R&L 1970: Lathrop
2000's dig "fixes a quadrature that does not correctly integrate Legendre
polynomials" is R&L's own closing sentence, "<3 % … S₈ Gauss weights". M&M 1984:
their entire accuracy concession — *"our scheme is only first-order accurate"* —
is one clause on the **second page**, in the sentence that contrasts them with
R&L.)

⭐ **Sharpening (2026-08-11) — the necessary-vs-sufficient case: COUNT THE
CONDITION'S RANK against the unknowns it is said to determine.** A secondary
routinely compresses a primary's derivation into "condition X determines Y", and
the compression **inverts the implication** when X is a *consequence* rather than
the *defining* requirement. Bailey-Morel-Chang: *"Forcing this β factor to be zero
determines the Morel and Montry weights."* β is **one scalar**; the weights are an
**N-vector** — the claim cannot be true, by dimension alone. M&M's real condition
is per-ordinate edge coincidence (N equations, their Eq. 15 + 16b) and β = 0 is a
parity corollary (their 17a–19). Cost of believing the gloss: you conclude any
τ with β = 0 is admissible and grant yourself licence to clamp/retune τ.
**The check is two mechanical steps**: (i) count equations vs unknowns in the
stated condition; (ii) if under-determined, *exhibit a second solution* — here,
root-finding τ₁ on Gauss-S8 with the other seven randomised gave β = O(1e-16) at
‖τ − τ_MM‖_∞ ≈ 0.24–0.31 (L-011's reproduce-it discipline, applied to a
uniqueness claim rather than a value). Then go find the primary's actual
determining condition; it is always stronger and always elsewhere. See
[[reed-lathrop-1970-angular-truncation]],
[[morel-montry-tau-angular-cell-edges]].

## L-014 — Two sources "conflict"? Check the CITATION EDGE first — the later one has almost certainly already ruled

`[M]` 2026-08-11. A brief posed a genuine-looking impasse: BMC Eq. 43 (τ barycentric
in μ) and Reed–Lathrop Eq. 16 (second order iff τ = ½ + O(w)) give **opposite**
answers at a level's grazing ordinates, "so we must choose, and we want to know
whether Morel & Montry saw this". They did — and the answer was not in a derivation.
**Morel & Montry cite Reed & Lathrop as their reference 1**, reject them in a
half-paragraph of the Introduction (R&L's induced quadratures cannot integrate
degree > 3, which breaks conservation under anisotropic scattering), and concede
their own order in the very next sentence: their scheme works with **any** quadrature
but is *"only first-order accurate"*. The whole "conflict" is a **declared trade**,
stated on the second page of the primary, decades before the secondary that appeared
to contradict it.

**Why:** an apparent conflict between two papers in the same lineage is nearly always
a *documented design choice*, because the later author read the earlier one and had
to justify diverging. Treating it as an open scientific question invites the agent to
adjudicate it itself — from numerics, from a third source, or worst, from taste — and
whatever it picks arrives with no provenance while a one-sentence primary ruling sits
unread. The failure is silent: an invented adjudication reads exactly like a sourced one.

**How to apply:** before adjudicating any two-source disagreement, spend one grep.
(a) **Open the later paper's reference list and look for the earlier paper.** If the
edge exists, (b) find the in-text citation — it is usually in the Introduction, in the
paragraph that says why this work is needed — and read the two sentences around it.
(c) Report the ruling as the *authors' trade*, naming what each side buys and pays,
not as your verdict. If the edge does NOT exist, say so explicitly: the conflict is
then genuinely open, and *that* is the finding. Companion to L-013 clause (c) — same
paragraph answers both questions.

## L-006 — Local literature folder FIRST, then Zotero, then Tier-2

`scratch/literature/` holds the user's full Nuclear Science &
Engineering run plus acquired PDFs. Check it (exact path) before any
online pull; the brief from the main agent should already name the
path. Zotero is Tier 1 after that; web databases are Tier 2.

**Why:** reading a local PDF is unambiguous and cheap vs re-acquiring;
pivoting to a *secondary* source is a structural decision needing user
approval, not agent autonomy (codified in `.claude/rules/delegation.md`).

**How to apply:** if a requested paper is not in the local folder,
report "not in local folder; acquire it, or will you add it?" — do NOT
unilaterally substitute a different source.

## L-008 — Mistral-OCR table BODIES live in the .mocr.json, not the sidecar

The sidecar markdown renders a table as a placeholder link
(`[tbl-0.md](tbl-0.md)`); the actual markdown table content is in the
raw cache at `pages[k]['tables'][j]['content']` of
`scratch/literature_ocr/<stem>.mocr.json`.

**Why:** grepping only the sidecar makes a tabulated report look
value-free (LA-3186's 21 tables were all placeholder links); render-only
transcription of dense number tables is slow and content-filter-risky.

**How to apply:** grep the sidecar to LOCATE tables (`tbl-`), then dump
bodies with a 5-line json loop over `pages[*].tables`; STILL verify
load-bearing values on the rendered page (OCR dropped/garbled cells in
LA-3186 Table I's n=20 rows). See [[la3186-level-symmetric-quadrature]]
for the worked pattern. Suggest (user-owned tool): `ocr_literature.py`
could inline `tables[].content` into the sidecar at emit time.

## L-009 — "403" is not "paywalled": check the OA STATUS before telling the user a paper is inaccessible

A `403` from a publisher is usually a **Cloudflare bot challenge**, which a
human browser walks straight through. Before reporting a paper as
inaccessible, resolve its OA status with Unpaywall
(`api.unpaywall.org/v2/<doi>?email=…`) or OpenAlex `open_access.oa_status`.

- **`bronze`** = free-to-read on the publisher's own site, no licence. The
  paper IS obtainable at zero cost — say "download it in a browser", NOT
  "paywalled". Its tell in OpenAlex is an `oa_url` carrying
  `?needAccess=true`: that reads like a false-positive OA flag and is not one.
- **`green`** = a repository copy exists; `locations[]` names it. Conversely
  `any_repository_has_fulltext: false` + a single `locations[]` entry means
  **no green copy exists** — stop hunting institutional repositories.

**Why:** the two verdicts trigger opposite user actions (one click vs an
acquisition decision), and the wrong one either wastes the user's money or
strands a free paper. Measured 2026-08-11 on Yamamoto 2007 JNST 44(2):
`403` on every T&F route, yet `bronze / publisher / publishedVersion`.

**Also check whether the archive MOVED.** J-STAGE 404s on all JNST
`article/jnst/44/*` paths because JNST vols. 1964-2011 were migrated to
Taylor & Francis; the original `10.3327/jnst.*` DOIs now redirect to
`10.1080/18811248.*`. A 404 at the expected native host means "relocated",
not "absent" — follow the DOI. See [[ty-polar-quadrature-moc]].

## L-011 — When an equation's SYMBOLS are ambiguous, REPRODUCE THE PAPER'S TABLE; the numbers decide what no amount of re-reading can

`[M]` 2026-08-11, BMC 2010 Eq. (41)/(75). The β sum is written with symbols
`μ_{m±1/2}` that *look* like the quadrature's cumulative-weight cell edges. Read that
way, β is identically zero for every scheme (Lathrop's Eq. 57 oddness argument proves
it), which flatly contradicts BMC's own Table I. Re-reading the page cannot settle it —
both readings are typographically identical. **Coding the candidate reading and
reproducing Table I did**: the scheme-IMPLIED edge march (Eq. 43 solved forward from
the τ under test) reproduced all 8 printed values to every digit
(`7.698004e-01`, `2.06E-01`, `−3.57E-03`, …), the cumulative-weight reading gave
round-off everywhere.

**Why:** a paper's numeric table is a **specification of its own equations** — the one
artefact that pins the meaning of every symbol simultaneously. Prose can be ambiguous
and a scan can be faithful to an ambiguity; a 4-row table cannot be satisfied by the
wrong reading. This is also the only way to detect that two papers using the same
Greek letter mean different objects.

**How to apply:** whenever an extraction turns on *which* quantity a symbol denotes —
and especially when two sources appear to contradict each other — look for a tabulated
result and reproduce it in `.venv/bin/python` before reporting. Table bodies are in the
`.mocr.json`, not the sidecar (L-008). Report the reproduction as `[M]` evidence for
the reading. Corollary: **also compute the diagnostic on the USER's own rule** — the
same pass found BMC's β is round-off for *every* convention on a symmetric
equispaced-ω march, so recommending it as an oracle would have shipped a blind gate.
See [[morel-montry-tau-angular-cell-edges]].

## L-010 — Scan-verification proves the OCR faithful, NOT the equation correct: the PRINT itself can be wrong

`[M]` 2026-08-11, Bailey-Morel-Chang 2010 NSE 165. The sidecar rendered Eq. (52)
as `μ_{m+1/2,n} = μ_{m+1/2,n} + w̄_{m,n}` — self-referential, hence impossible. The
obvious diagnosis is an OCR slip; the rendered page (PDF p. 10) shows the **same
thing**. It is a **typo in the published journal**. The correct RHS subscript is
`m−1/2`, proven two ways: the sphere's analogous Eq. (12) is printed correctly
(`μ_{m+1/2} = μ_{m−1/2} + w_m`, verified p. 5), and the recursion is otherwise
vacuous.

**Why:** the `[VERIFIED-ON-SCAN]` marker answers "does the sidecar match the page?"
It does NOT answer "is the page right?". An agent that stops at scan-match will
transcribe a published typo into a solver with the strongest available provenance
badge attached — an ERR-032-class hazard wearing the ERR-032 countermeasure.

**How to apply:** for any load-bearing equation, run a **third** check beyond
OCR-vs-scan — cheap and mechanical, pick whichever applies:
(a) **the analogue** — curvilinear papers derive sphere and cylinder in parallel;
    cross-read the twin equation, since a typo rarely hits both;
(b) **non-vacuity** — a recursion whose LHS and RHS name the same symbol, a
    definition that defines nothing, a sum that cannot close;
(c) **the seed/terminus** — BMC state `μ_{M+1/2,n}` "will always be" `+√(1−ξ_n²)`
    *when computed recursively*; a mis-transcribed recursion breaks that.
Report the slip explicitly and give the corrected form with its justification.
Same pass also caught a genuine OCR fraction slip (Lathrop Eq. 30, `−⅓` in the
sidecar vs `−¼` on the page) — the two failure modes coexist and are distinguished
only by opening the page. See [[morel-montry-tau-angular-cell-edges]].

## L-012 — A "no source does X" finding needs its DENOMINATOR, and the named PRIMARY is usually the one you haven't read

`[M]` 2026-08-11. A prior deliverable of mine concluded *"the geometric arc-half-angle
edge, **which no source uses**"* — a negative existential. Re-audited: it was checked
against **2** sources, both SECONDARY. Widening to 3 found no counterexample, so the
claim survived — but the text Hébert himself uses to introduce the construction names
**Alcouffe & O'Dell** as its author (his ref. [36]), and that paper is **not local, and
4 databases could not resolve it**. So the authority for the very construction being
adjudicated had never been read, and the sentence gave no hint of that.

**Why:** an absence-of-evidence claim reads exactly like a positive finding, and it is
the one claim type whose confidence *should* scale with the denominator. Worse, the
sentence was ALSO imprecise: BMC's radial-cosine edges (Eq. 52) DO coincide with the
geometric arc values whenever the quadrature makes weight = cell measure — so
"no source uses geometric edges" conflated two different face quantities and hid a
quadrature-dependent equivalence.

**How to apply:** three mechanical steps before shipping any "no source …" line.
(a) **State the denominator** — "3 of 3 local sources that define a face cosine",
never a bare "no source". (b) **Grep the local text for who it CREDITS** — a chapter
that says "the approach proposed by X" is telling you the primary is elsewhere;
resolve or explicitly flag it as unread. (c) **Restate the claim as a MECHANISM, not a
survey** — "all three define the face cosine by a conservation recursion; whether that
lands on the geometric value is a property of the QUADRATURE" is checkable and
transferable, where "nobody uses X" is neither. See
[[morel-montry-tau-angular-cell-edges]].

⭐ **Sharpening (2026-08-12) — the denominator can be silently ZERO because the
INSTRUMENT is blind to the domain. Before trusting a search's emptiness, verify the
tool returns TRUE POSITIVES on the topic.** `[M]` Q68: I wrote *"the (η,φ) 2-D angular
treatment — searched hard, found NOTHING"* on the strength of a corpus grep plus
**seven** OpenAlex full-text formulations. The OpenAlex runs were worthless: **`SN`
matches *supernova***, so every query returned WMAP, Type Ia spectropolarimetry, MHD
and sea-ice papers. Zero of ~70 hits were transport papers. Re-running **one**
ISSN-scoped CrossRef `journal_search` on NSE found **Chaland & Samba 2016** — a full
2-D `(μ,φ)` product-mesh angular closure — on the first query, and the paper was in the
user's own archive all along.

**Why:** a search that returns junk and a search that returns nothing *look identical
in the negative*, because both end with "no relevant hits". The instrument's failure
is invisible in its output; it is only visible in the hits you DIDN'T read. And a
false negative is the most expensive error class here — it closes a question and tells
the user to go do original research.

**How to apply, three cheap steps.** (i) **Sanity-probe the instrument**: query it for
a paper you KNOW exists in the target domain; if that misses, the tool cannot see the
domain and its silence is meaningless. (ii) **Watch for homograph collisions in the
key token** — `SN`, `CP`, `MC`, `P1`, `DG` are all overloaded across physics; prefer a
**venue-scoped** query (CrossRef `journal_search` by ISSN) over free text whenever the
subject is nailed to a small set of journals. (iii) **Prefer the CITATION GRAPH for
"what came after X"** — Semantic Scholar `get_citations` on a seed DOI found every
post-2010 follow-up that title search missed, and its denominator is *countable and
reportable* ("Lathrop 2000 has 4 citers, all classified"), which is exactly what
clause (a) demands. Companion: [[user-nse-volume-archive]] (the *other* half of the
same failure — the local denominator was also bigger than I assumed).

⭐⭐ **Second worked case (2026-08-25), and the cost profile that makes this the
worst lesson to repeat: MY OWN "#158 curvilinear LD is unpublished" negative had
shipped into PRODUCTION docstrings** (`linear_discontinuous.py`, `scheme.py`) as a
guard message + capability-flag rationale, telling every future session to derive
from scratch what Adams-Martin 1992 NSE 111 App. A prints in full — and that paper
was IN `scratch/literature/` when re-searched. Two mechanics worth keeping:
(a) **a refuted negative does not stay in the memo that made it — it propagates
into code and issues; the correction owes a sweep of those surfaces** (flag them to
the main agent; a literature agent does not edit production). (b) **Re-run the
corpus grep when the corpus GROWS**: the refuting PDF's sidecar was a later-added
`.textlayer.md` variant (Mistral 401'd 2026-07-26, different suffix), so any
earlier grep with `--include='*.md'` narrowed to mocr stems, or run before the
add, missed it silently. A negative dated before the corpus's newest file is
STALE, not settled. Full record: [[issue-158-linear-discontinuous-cell-update]].

## L-007 — Recognize a dead Zotero server and fail over to Tier 2 immediately

A Zotero MCP server that returns 0 hits on known-present items together
with connection-refused on port 23119 is BROKEN, not empty.

**Why:** repeated 0-hit Zotero calls waste turns; the pattern is
diagnostic. Several recent sessions ran with Zotero DOWN.

**How to apply:** on the 0-hit + conn-refused signature, stop querying
Zotero, note "Zotero down this session — no annotations checked" in
output, and use the local folder + Tier-2 databases. See
[[reference-zotero-flakiness]].

## L-015 — A SPEC question is not a literature hunt: the primary standards are one `curl` away, and `grep` cannot read what you extract from them

Two failure halves, both measured 2026-08-31 while establishing the ENDF-6 / GENDF
(n,2n) format rules.

**(a) Don't Tier-2 a standard.** `scratch/literature/` is a *transport-theory* library —
78 PDFs, **zero** nuclear-data-format documents, and a wide-vocabulary grep
(`ENDF|NJOY|GROUPR|GENDF|Kalbach|LAW=|ACER`) returned only 6 incidental mentions inside
physics papers. But the **ENDF-6 Formats Manual** (CSEWG ENDF-102 / BNL-203218-2018-INRE,
418 pp) and the **NJOY2016 manual** (LA-UR-17-20093, 816 pp) are freely published,
**born-digital with a real text layer**, and download in seconds from fixed URLs. Getting
them IS getting the primary source. A whole class of question — *"what does the FORMAT
permit?"* — is answered by two `curl`s, not by a database search. URLs, page offsets and
the settled rulings: [[endf6-gendf-njoy-n2n-formats]].

**(b) ⛔ `grep` here is ugrep with `-I`, and it SILENTLY SKIPS a pypdf-extracted `.txt` as
binary — 0 hits, no error, identical to a clean file.** I lost a cycle believing the
extraction had failed. This is the L-family filter defect in a new dress (a broken filter
and an empty corpus print the same thing), and it is worse than the `\bword\b` case
because nothing about the pattern looks wrong. ⟹ **extract with page markers and search in
PYTHON**, never the shell:
`parts = re.split(r'@@@PDFPAGE (\d+)@@@', txt); pages = {int(parts[i]): parts[i+1] for i in range(1,len(parts),2)}`.

⚠ **Rider, and it upgrades L-010 rather than repeating it: `pypdf` breaks FRACTIONS across
lines.** `(2ℓ+1)/4π` extracted as `2ℓ + 1 / 4π`. On a normalization question that is the
*entire* answer — `/2` and `/4π` differ by the `2π` azimuthal integral. So the
rendered-page check is mandatory **even for born-digital PDFs where L-010's OCR hazard does
not apply**: the new mechanism is *layout loss*, not character error. Both load-bearing
equations were verified visually and both held.

## L-016 — When the manual under-determines the behaviour, MEASURE THE DATA — a census of real files outranks an unprinted table

`[M]` 2026-08-31. The question was *"does GENDF's `NL` reflect what the evaluation
contains, or only what the user requested?"* The NJOY manual gets you close and then
**stops**: it says `lord` is a ceiling, names `subroutine init` as choosing a per-reaction
default, and **does not print `init`'s table**. Honest answers available from the manual
alone: "≤ lord+1, mechanism unstated".

What settled it was **fetching the nine actual ENDF/B-VIII.0 evaluations behind the user's
nine data files and parsing how each represents the reaction** (~40 lines of Python). The
answer fell straight out and was *stronger* than the manual could have been, because the
corpus contained its own controls: one nuclide's file is a **two-line declaration of
isotropy** and gets `NL=1` under the *same* processing run that gives every other nuclide
`NL=7`. A counterexample inside the user's own corpus is unimpeachable evidence.

**How to apply.** When a brief asks "is my measurement a property of the format, the
evaluation, or the processing?", that is **three questions with three different sources**:
the spec (what is permitted), the evaluated file (what this nuclide does), and the
processor's manual + the processed file (what was done to it). Get all three; they are
usually all obtainable, and the disagreements between them are the finding. Corollary
worth its own line: **reproduce the user's headline number from the raw upstream data.**
Computing µ̄ independently off the ENDF tape (+0.25…+0.27) and off the GENDF bytes
(+0.2783) turned "your measurement is consistent with the spec" into "your measurement is
confirmed at three independent levels" — and it cost one script. Companion to L-011
(reproduce the paper's TABLE), one tier down: reproduce the user's MEASUREMENT.

⚠ And the caveat that made the deliverable useful rather than merely affirming: the same
census showed **8 of the 9 nuclides are `LCT=2` (centre-of-mass)**, where the processor's
CM→LAB transform manufactures ℓ≥1 moments *from an isotropic distribution* — so the
confirmation transfers to exactly ONE of the user's nuclides. **A confirmation owes its
denominator just as much as a negative does** (L-012's clause, applied to a positive).
