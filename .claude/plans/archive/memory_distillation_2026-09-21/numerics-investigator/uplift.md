# Uplift candidates — numerics-investigator, 2026-09-21

Eleven meta-lessons that generalise beyond this agent. Each names the rule or skill it would
join, proposes a clause in that artefact's own form (a `check:` and a `tell:` where the artefact
uses them), and names the founding case. **I propose the clause; the orchestrator edits the
source page** (`docs/development/rules/*.md` or `docs/development/skills/*.md` — the
`.claude/` copies are generated).

Ordered by leverage. U1–U4 are the ones I would take first.

---

## U1 — `vv-principles` anti-pattern #6 gains the project's-own-corpus tell

**Why.** #6 ("NEVER trust a reference not traced to a structurally-independent analytical/
symbolic ground") lists three tells, and all three are *solvers*: "MC vs MC; CP vs unverified
MC; method-of-images converged to the wrong BC". Nothing warns that **a theory page and the
docstring that quotes it are ONE source, not two** — which is the shape that actually bit.

**Proposed addition to #6.**
> **check:** when the reference is a project artefact, ask what its author read. A theory page
> and the docstring quoting its equation are ONE source; so are a page and the test whose
> expected value was copied from it. Trace to the literature equation number or to a sibling
> module of a DIFFERENT kernel family.
> **tell:** "the code matches the theory page" offered as verification.

**Founding case.** ERR-063: `docs/theory/peierls_nystrom.rst` Eq. `peierls-mg-operator` AND
`solve_peierls_mg`'s docstring math both documented the same wrong `χ_g(r_i)` sink index.
Digest L8.

---

## U2 — `numerical-bug-signatures` Signature 8 gains the vacuous-tolerance-sweep discriminator

**Why.** Signature 8's first discriminator is "Tighten `inner_tol`: the fixed point does NOT
move ⟹ it IS a discretization bug". That inference is **false whenever every swept tolerance is
tighter than the residual the capped run reaches** — then all runs hit `max_iter`, return the
same iterate, and tolerance-insensitivity reads as a floor when it is the cap. The signature
currently walks its reader into the wrong verdict.

**Proposed addition to Signature 8's "Discriminator" block.**
> ⚠ **A tolerance sweep is informative only if at least one swept tol is LOOSER than the
> residual the capped run actually reaches.** Read `history.converged` and `n_inner` against
> `max_inner` FIRST: a plateau at exactly `max_iter − 1` is the whole diagnosis, and
> bit-identical answers across four decades of `inner_tol` are then the CAP speaking, not the
> discretization.

**Founding case.** A gate red at 3.29e-10, bit-identical across `inner_tol ∈ {1e-9 … 1e-15}`
while the running residual at the cap was 1.185e-09. Digest L10b.

---

## U3 — `vv-principles` anti-pattern #5 (or the bit-identity section) gains "measure the order at the cell"

**Why.** #5 says "convergence rate is not convergence value". The dual failure has no clause: a
**localized** defect measured in a **volume-weighted** norm is diluted below the floor
(`√V ~ h^1.5` at a pole cell), so the global L2 reports a clean order while the cell is O(h).
The rule's own bit-identity section demands an end-to-end swap, a silent control and an AMPLIFY
leg — and all three can pass while the mechanism stays invisible in the norm.

**Proposed addition.**
> **check:** for a defect localized to a boundary row, a pole cell or one interface, measure the
> convergence ORDER at that cell, not in the norm — a volume- or measure-weighted norm dilutes
> it by the cell's own weight. If the cell is already at the target order with the suspect term
> disabled, the term cannot be repairing a deficiency that is not there.
> **tell:** a global L2 order table offered as evidence about a single-cell closure.

**Founding case.** The curvilinear τ-clamp / pole-floor thread and the LD boundary-slope
verdict. Digest L12.

---

## U4 — `numerical-bug-signatures` gains a new Signature: the greedy-`Ellipsis` spectator axis

**Why.** It is a recurrent AI-authored array-indexing bug class with a sharp 2-D fingerprint,
and the skill's own "Add-a-signature protocol" is satisfied: the class recurred (four verbs in
one commit), it has an exact diagnostic probe, and it names the test classes blind to it.

**Proposed signature (abridged to the skill's template).**
> **Symptom:** a moment-tensor / einsum path is bit-identical for the single-moment (scalar)
> case but raises `IndexError` on a rectangular grid, or returns a SILENT wrong value on a
> square grid with an asymmetric material map, as soon as a trailing spectator axis (LD's `2^d`
> spatial-moment axis) is present.
> **Mechanism:** a per-cell fancy index `cells = (Ellipsis, *idx)`. `Ellipsis` is greedy from
> the front: it absorbs the leading `(m, g, …)` axes, so `*idx` lands on the LAST k axes, which
> under a trailing axis are `(…, last-spatial, trailing)` instead of the two spatial axes.
> **Diagnostic probe:** run the same path with `spatial_moments=2` on a RECTANGULAR (`nx ≠ ny`)
> flux; isolate by monkeypatching the verb back to `(Ellipsis, *idx)`.
> **Blind test classes:** every scalar-moment test — with no trailing axis
> `Ellipsis ≡ (slice, slice)` and the two spellings are bit-identical.
> **Fix:** pin the leading axes explicitly, `cells = (slice(None), slice(None), *idx)`.
> **Failure mode:** #2 (variable/axis swap), gated by a spectator-axis presence.

**Founding case.** #276 A2, commit `0b3275d`, all four `MaterialXSField` moment-scatter verbs.
Digest L13. (No ERR entry exists — the protocol asks for one first, so the orchestrator may
prefer to file the ERR and then the signature.)

---

## U5 — `vv-testing` gains a fixture-consistency clause

**Why.** `vv-testing` governs tolerances, markers and levels but says nothing about the fixture
DATA being internally consistent. An inconsistent mixture makes two legitimate references
disagree with **no bug in either**, and the investigation that follows is unfalsifiable — I
spent a probe chasing a "benign pole" that did not exist.

**Proposed clause, under a new heading "A hand-built mixture is gated on its own identity".**
> A hand-built verification fixture prints and asserts the mixture's own consistency identity
> `σ_t == σ_c + σ_f + Σ_to SigS[0][g,:]` before any reference value is trusted.
> **check:** an inconsistent mixture gives the transport balance (removal `σ_t`) and the
> production/absorption balance two DIFFERENT answers, and two correct solvers will report
> different ones; a `[to, from]`-vs-`[from, to]` transposition of `sig_s` is the common cause.
> Check the companion poisons in the same line: a group with `φ ≡ 0` is a 1-group problem
> wearing a 2-group shape (anti-pattern #3), and a negative `σ_c` is unphysical.
> **tell:** a brief's reference value on a hand-built mixture, quoted but never re-derived.

**Founding case.** #340 N5: the brief's "benign pole" was 30 % off; one character (`sig_s=s.T`)
repaired it and the intended pole then appeared. Digest L20·5.

---

## U6 — `instrument-doctrine` X1 gains the TRANSFER GAIN as a design-time procedure

**Why.** X1 asks "name the input that would make this read differently". The prior question — *can
this statistic bound that quantity AT ALL?* — is answered by one measurement that costs less
than any threshold hunt, and the skill's procedures do not name it. Every "let us add a
certificate / a health metric / a tolerance on X" proposal is this question.

**Proposed addition to the skill's X1 section.**
> **Transfer gain, before any threshold.** For "can statistic `X` gate quantity `y`?", measure
> `|Δy| / X` across configurations FIRST. A threshold on `X` bounds `|Δy|` only through that
> gain, so an unbounded gain means no constant exists and the tuning exercise is void. Report
> the gain's SPREAD and the two populations' overlap, not a candidate threshold.
> **tell:** a proposed threshold with a sensitivity/false-alarm table and no gain column.

**Founding case.** #340 N5: `|Δk|/defect` spanned 1.16e+05×, populations overlapped 634×, and a
zero-false-alarm threshold missed 15 of 16 corrupting cases. Digest L20·1.

---

## U7 — `vv-principles` bit-identity section gains "report the relative error before the nulp"

**Why.** The section already prescribes `assert_array_almost_equal_nulp(nulp=K)` and bounds the
legitimate drift as `reduction depth × ULP`. It does not warn that **the nulp reading itself is
uninterpretable in both directions** — the received wisdom covers only "huge nulp near zero is
nothing", and the dual (a huge nulp that means a percent-level difference) is what actually
mis-sorts a triage.

**Proposed addition.**
> **check:** on any bit-identity red, print `max|a−b| / max|b|` BEFORE reading the nulp count.
> `[M]` `1.04e+15` nulp has meant an 8 % difference over 216/216 elements, and an
> `array_equal → False` on two arrays that print identically has meant 1 ULP (`2.06e-16`). The
> relative error sorts the triage; the nulp count cannot.
> **tell:** a nulp figure quoted as the magnitude of a regression.

**Founding case.** The 9-red stale-reference triage. Digest L22·1.

---

## U8 — `vv-principles` gains the tolerance-pin / bit-identity-pin gap

**Why.** A correct rule-tier tolerance pin and a correct consumer-tier bit-identity pin can be
simultaneously satisfied and simultaneously wrong about a change: a 3-ULP move is INSIDE the
producer's `< 8 ulp` contract and OUTSIDE every downstream `array_equal` / `sha256`. The gap is
structural, not an oversight, and nothing in the corpus names it.

**Proposed clause (a new anti-pattern, or a paragraph under bit-identity).**
> **NEVER** expect a producer's tolerance pin to warn about its consumers' bit-identity pins —
> **instead** land a byte-level fingerprint beside the tolerance pin whose sole job is "the
> bytes moved, re-baseline the consumers".
> **check:** for every primitive gated at a tolerance, grep its downstream `array_equal` /
> `assert_array_equal` / digest consumers; if any exist, the producer owes the fingerprint.
> **tell:** a producer gated "correctly, since neither construction is *the* answer" with
> frozen-byte consumers downstream.

**Founding case.** `gauss_legendre` vs numpy's `leggauss` at `< 8 ulp`, five downstream
snapshot rows red. Digest L22·3.

---

## U9 — `code-search` gains two probe-harness silent-zero mechanisms

**Why.** `code-search` already carries two silent-zero cases (the ugrep anchor-in-alternation,
and `git grep -E` having no `\b`) and the rule that a completeness claim needs a positive
control. These two are the same family, same failure shape (a confident, empty, wrong answer),
and both bit me in one session while measuring across a worktree.

**Proposed additions, in the page's existing form.**
> - **`grep` against pytest output reads through ANSI colour codes.** `grep -cE "^FAILED"`
>   reported `0` at nine commits including two already measured RED, because pytest's colour
>   escapes precede `FAILED`. check: pass `--color=no` to pytest, or drop the `^` anchor.
> - **`python -c` prepends the CWD to `sys.path[0]`, ahead of `PYTHONPATH`.** A worktree probe
>   run as `python -c` from the main tree silently imports the MAIN tree and prints HEAD's
>   values for every commit under test. check: run probe SCRIPT FILES located outside the
>   repository, and make every probe print `module.__file__` as its first line.
> - tell: a per-commit sweep whose numbers are suspiciously identical across commits.

**Founding case.** The 9-commit bisect in the stale-reference triage. Digest L22·5.

---

## U10 — `instrument-doctrine` X2 gains "parameter-independence is a derivation hint"

**Why.** X2 governs populations and denominators; this is its constructive twin. When a fitted
law is measured INDEPENDENT of a parameter the object plainly contains, the fit is the wrong
instrument — the governing equation is combinatorial and a closed form is reachable, usually in
minutes where the numerical route is hours.

**Proposed addition (X2's census protocol, or a short new paragraph).**
> **check:** before shipping a fitted law, list the parameters the object contains and the ones
> the law does NOT depend on. Each independence is a constraint the closed form must satisfy and
> a hint that the governing equation is combinatorial: substitute the degenerate branch of the
> object's own closure and see what cancels.
> **tell:** a counting law with a fitted coefficient and no derivation, whose measured
> independence list is longer than its dependence list.

**Founding case.** #344: `dim ker` measured mesh-, `c`- and cross-section-independent and
exactly `∝ ng`; the closed-form basis then fell out in `0.05 s` against a `23 s` dense SVD at
half the size. Digest L24·1.

---

## U11 — `vv-principles` #24 gains a clause (f): a threshold has two edges

**Why.** #24(d) is the zero-set check ("solve `instrument = 0`; if the solution IS the incumbent
the instrument measures distance-to-incumbent"). The adjacent failure is that the justifying
instrument is **MONOTONE** in the threshold, so a flat scan is a true reading that carries zero
information about one of the two edges — and the guard that forecloses the other edge is usually
already in the code, undocumented.

**Proposed clause (f).**
> **(f) *two edges*** — a numerical threshold bounds a quantity from both sides, and the
> instrument offered as its justification is usually monotone in the threshold, hence
> structurally blind to one edge. Before citing a flat scan, ask which direction the statistic
> can even move in; then find the guard that already forecloses the other side and BISECT its
> refusal boundary rather than reasoning about it.

**Founding case.** `_DENSE_METRIC_RCOND` (ERR-080 / #429): the Parseval ratio read a flat
`1.000000000` across `[1e-15, 1e-2]` while a construction-time pair-consistency guard refused
everything below `8.696754e-17`. Digest L27·2–3.
