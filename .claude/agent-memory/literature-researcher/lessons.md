# Literature Researcher — Lessons

Behavioral corrections only (failure → rule). The AGENT.md / `research`
skill hold source priorities, notation mapping, and extraction
procedure — never duplicate those here. Durable paper-extraction
RECORDS (equations, notation maps, value tables) live in the topic
files indexed by `MEMORY.md`, NOT here.

---

## L-001 — Bailey is the canonical cylindrical SN reference; state the geometry when mapping µ

Bailey, Morel & Chang (2009) NSE: **Eq. 50** = α recursion, **Eq. 74**
= Morel-Montry weights — the exact equations ORPHEUS implements for
curvilinear discrete ordinates.

**Why:** Bailey's µ = radial cosine = ORPHEUS `mu_x`, but in
slab/spherical contexts µ means the *streaming* cosine. The symbol is
overloaded across geometries.

**How to apply:** ALWAYS name the geometry when mapping µ. For
curvilinear SN balance / the ΔA/w geometry factor, start at Bailey §4,
not Lewis & Miller (who omit the asymptotic analysis explaining *why*
the geometry factor exists).

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

## L-007 — Recognize a dead Zotero server and fail over to Tier 2 immediately

A Zotero MCP server that returns 0 hits on known-present items together
with connection-refused on port 23119 is BROKEN, not empty.

**Why:** repeated 0-hit Zotero calls waste turns; the pattern is
diagnostic. Several recent sessions ran with Zotero DOWN.

**How to apply:** on the 0-hit + conn-refused signature, stop querying
Zotero, note "Zotero down this session — no annotations checked" in
output, and use the local folder + Tier-2 databases. See
[[reference-zotero-flakiness]].
