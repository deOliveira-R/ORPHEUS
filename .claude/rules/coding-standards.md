# Coding standards — the minimum-quality floor

These are **preemptive minimum standards** every contributor (main agent and every
sub-agent) follows by default, regardless of whether they are chasing elegance. They are
the *floor*. For the *ceiling* — recognizing how bad code manifests and refactoring toward
excellence — code-producing agents load the **`coding-elegance` skill**. (See also Cardinal
Rules 1 "Correctness" and 2 "Architecture" in CLAUDE.md.)

## Clean before extending

Before adding a capability to a class/module — especially on a new design — first run a
cleanup pass on that layer: collapse double paths, move concepts to their native place,
delete dead shims, fix twin sources of truth. The new capability should then land as a
**no-op extension through the single generic body**, not as a new arm grafted onto
structural debt.

- When a plan proposes a capability extension, insert a dedicated **cleanup phase before
  the capability phase**. Order findings into must-precede vs independent-polish vs
  explicit-WAIT; gate each cleanup substep as bit-identical where possible.
- Rationale: extending a layer that carries debt grows a third arm on every seam. (C5,
  2026-06-11: `from_axes` round-tripped axes→legacy-mesh→axes, so 3-D admission would have
  needed a new arm in the converter, constructor, AND trace gate — until the cleanup
  inversion let 3-D flow through the one generic body.)

## Type vs property — before minting a type

A representation earns its own **type** only when it is genuinely a different object, not the
same object wearing a label. The decidable (grep-checkable) criterion: mint a type **iff**
(a) there are **≥2 non-isomorphic bases/realizations** of the concept, AND (b) a
**non-identity morphism** is actually applied to it. Otherwise the concept is a **property**
— a field or flag on an existing type.

- If the only "change of basis" is the identity (one realization, no transform), a separate
  "type" is theatrics: it adds ceremony and a conversion seam without making any illegal
  state unrepresentable. Make it a property.
- This is the *type-minting* corollary of the `coding-elegance` "defer abstraction until ≥2
  instances" rule. Worked: an expanded-order spatial moment with a single basis and an
  identity change-of-basis is a `property`, not a `SpatialOrder` type.
- **Corollary — an axis that changes the ARITHMETIC INTERFACE cannot be a phantom type
  parameter.** The criterion above decides *whether* to mint; this decides *how to encode
  it*. `Generic[Tag]` is erased at runtime and does not specialize dunders, so every
  instantiation shares ONE `__add__` body. If two values of the axis need different
  arithmetic signatures — a torsor `A×V→A` that must FORBID `A×A`, versus a vector
  `V×V→V` — no shared body satisfies both, and the encoding must be a distinct class.
  Decision lattice: axis changes the **arithmetic** ⇒ class; axis changes the **shape** ⇒
  class; only an axis that changes NEITHER may be a phantom parameter. Negative
  discriminator: an implementation that "passes" only by branching on a stored tag field
  at runtime is stringly-typed dispatch — `replace(obj, tag=Other)` type-checks and walks
  straight through the gate the type was minted to be.

## Retire as you go (aggressive retirement)

Superseded code is **noise that obscures signal** — it makes the codebase harder to read,
breeds "which path is canonical?" confusion, and invites accidental extension of the wrong
path. Retirement is a first-class deliverable, not optional cleanup.

- Every refactor that introduces a more elegant pattern **MUST retire its predecessor**.
- Deprecation aliases / compatibility shims live for **one merge cycle only** — remove them
  on the next pass. Never preserve backward-compat code unless the user explicitly
  authorizes it (e.g. a public API on a versioned release).
- Treat the **retirement audit as its own numbered substep**; a plan should carry a
  "retirement list" enumerating what gets deleted (with `file:line`).
- **A retirement's own past-tense NOTE is a confidence trap — discriminate the surviving
  references BY TENSE.** A batch that deletes a symbol and adds an honest "`X` existed
  until 2026-07 and was retired with zero consumers" note *reads* like a completed audit;
  it is not. Grep the deleted name across the whole tree and sort the hits: past-tense
  ("existed", "was retired") is correct history and STAYS; a **present-tense claim**
  ("provides an ergonomic shortcut") or an **imperative instruction** ("Apply the marker —
  `@verify.lN(...)`") is a MUST-FIX — a contributor following it hits `ImportError`, and a
  maintainer re-adds the symbol "to match the doc", reopening the retired twin. Expect the
  offenders to be PRE-EXISTING lines outside the diff, sometimes ~50 lines from the
  batch's own note.
- **Retirement means test migration:** retiring a symbol includes **rewiring its tests to
  the successor**, not deleting them with the symbol. Behavioral test (correctness contract)
  → rewire to the new API; API-smoke test (symbol exists) → delete; characterization test
  (e.g. FD-vs-WDD delta) → keep under `tests/<module>/characterization/`. Pure delete-only
  retirement is incomplete — it loses coverage. Inventory with `grep -rn "<symbol>" tests/`.
- **Retirement means MARKER migration too.** A retired test takes its
  `@pytest.mark.catches(...)` / `verifies(...)` with it: the successor asserting the same
  invariant must be re-tagged, or the coverage edge silently disappears while the
  error-catalog entry keeps naming the dead test class — an audit "MISSING" whose stated
  L0 test still *reads* plausible. Grep the catalog and the `tests/_harness` registry for
  the retired symbol alongside the code grep.
- **A rewire can silently DEMOTE a gate's claim class without touching one line of the
  test body.** When a retirement re-points a comparison target at the successor, re-ask
  **"are the two sides still INDEPENDENTLY produced?"** If the survivor is the *caller* of
  the other, a two-implementation bit-identity gate has become a value compared with
  itself through a wrapper: green forever, keeping its authoritative name, unable to
  detect the drift its docstring advertises — and invisible in review because BOTH sides
  are genuine production calls. The tell in a diff is a local variable still called
  `legacy` beside a brand-new API. The mutation test is one line: replace the SUT's body
  with garbage; if the pin stays green it was never a pin. Do not delete the gate —
  re-scope every doc and docstring crediting it, and name the pin that actually survives.
  **Then check that the named replacement is real:** a redirect to "the regression
  snapshots pin this" is worthless if those snapshots were re-baselined BY the same carve
  (measured 2026-08-03: three `cyl_*` snapshots cited as the pre-carve anchor had all been
  re-captured by the consolidation commit itself).
- **The retirement audit's blast radius is THREE searches, not one** (4–5 agents converged on
  this independently): (1) **graph callers** (`nexus impact`/`callers`) — necessary but NOT
  sufficient; the call graph misses property-reached leaves (`callers()==0` but live via a
  `cached_property`), class-name *bypass* consumers, and direct constructors of a guarded
  type; (2) **text-grep the symbol across code, tests, AND `docs/`** — an unresolved
  **Python-domain** cross-reference (`:func:`/`:class:`/`:meth:`/`:mod:`) renders as plain
  text with **no `-W` warning**, so the Sphinx gate does NOT catch a code retirement's doc
  blast radius — **and `-n` (nitpicky) does NOT save you either.** (This clause read
  "unless the build runs `-n`" until 2026-08-03; that was false, and it told every
  retirement audit it was covered when it was not.) Sphinx can only nitpick what it
  RENDERS, so a docstring in an un-`automodule`'d module is invisible at EVERY severity,
  as is every file under `tests/`. That is the majority case here, not an edge case:
  the doc source carries only ~45 live `automodule` directives, and the whole of
  `numerics/measure.py`, `numerics/operator.py` and `numerics/quadrature/` is among the
  many with zero — which is exactly why a module retirement left 22 dead `:class:`/`:mod:`
  refs that no build of any severity could see. Before concluding "`-n` would have caught
  this", check whether the module is rendered at all; if it is not, **grep is the only
  gate**, and an unchanged warning count proves nothing about it. (Measured 2026-07-15, Sphinx 9.1.0:
  `:doc:` and `:ref:` **do** warn — `ref.doc` / `ref.ref` — so *page* moves and *label*
  retirements ARE gated by `-W`; the silent class is the Python-domain roles, plus **raw path
  strings** in prose/docstrings, which no build ever checks. A path assembled from segments —
  `REPO_ROOT / "docs" / "theory"` — is invisible to a path-grep too; grep the **last
  segment**.) **Grep the CONCEPT, not only the symbol:** a field/flag is documented in two
  registers — by NAME, which greps, and by PARAPHRASE, which does not. A `list-table`
  column headed "Sweep-cycle flag" carries per-law values with no symbol in any cell; one
  audit's 7 exact hits missed 17 further cells. After the symbol grep, grep the
  hyphen/space variants of the concept the symbol names. And the paragraph that JUSTIFIED
  the retired thing inherits its wrongness — re-verify that prose against the replacement
  rather than only deleting the dead name from it. (3) **direct constructors** of any guarded type (a guard-at-source change
  reaches every `T(...)` caller, not just the factory path). Run all three, then retire.
- **Retiring a MESSAGE STRING: grep the SHORTEST distinctive fragment, never the full
  sentence.** An exception/log message is an API the moment a test pins it, and tests pin
  **substrings**. A grep for your own longer wording is strictly LESS sensitive than the
  consumer's pattern, so it returns a confident, empty, wrong answer. (2026-08-06, G6.3 step
  8.0: retiring `OperatorSum`'s inline check onto a shared helper reworded
  `"OperatorSum requires equal domains"`. The audit grepped `requires equal domains` — which
  matches the production line — and reported only the definition site. Two gates matched on
  `"equal domains"` alone and went red in the wide run; a third reference was prose. The
  correct pattern was the two-word fragment.) Corollary: when the audit does find pins,
  prefer **keeping the established vocabulary** over renaming it — the phrase is load-bearing
  provenance ("this guard fired, not some incidental raise"), and it is greppable precisely
  because it has not drifted.
- **Routing a call site through the ALGEBRA silently raises its operand requirement from
  "has the verb" to "IS the type" — and only duck-typed test doubles notice.** Re-spelling
  `f.apply(g.apply(x))` as `(f @ g).apply(x)` is arithmetic-neutral by construction, so it
  reads as a pure re-spelling; it is not. `@` needs `__matmul__`, i.e. a real operator,
  where `.apply` needed only an attribute. Production is usually unaffected (its objects
  come from a factory that already returns the type), which is exactly why the breakage
  surfaces in a *test* and looks like a broken test rather than a contract change. Fix it
  by making the surrogate honour the contract it stands in for; do NOT add a runtime guard
  for a case the type system now covers (that is the harmful-stub anti-pattern). (2026-08-06,
  G6.3 step 8: a `_NoTransposeLaw` stub with only `apply` hit
  `TypeError: unsupported operand type(s) for @`.)
- **And a retirement onto a SHARED helper moves the raise's provenance one frame out.** Any
  gate asserting the innermost traceback frame is now asserting the helper, which is
  reachable on behalf of *every* caller — so the gate silently widens from "this composite
  refused" to "something refused". Re-point it at the helper AND assert the CALLER frame (or
  an owner tag in the message), or the retirement demoted a provenance pin while leaving its
  name intact — the same defect class as the fuller-view/bit-identity demotion above.
- **Mass-deletes are retirements too — and untracked shadow-copies mask the breakage.** A
  "chore: mass-delete old diagnostics" sweep owes each file the same 3-search audit, with two
  sharpenings: (a) grep the **module/script NAME**, not only its symbols — a subprocess-worker
  import (`from diag_x import f` inside a `textwrap.dedent` worker string) is text the call
  graph never sees; (b) a diagnostic consumed by a tracked test is **production
  infrastructure** (the instrument behind a pinned baseline), never "shipped/falsified"
  debris. If a consumer stays green after the delete, check WHERE its import resolves: an
  untracked shadow (a scratch/ working copy, a stale `__pycache__`) can serve the import and
  keep the breakage silent until the shadow evaporates. (2026-07-13: `15486f66` mass-deleted
  `diag_cin_aware_split_basis_keff` while the CP rank-n protocol test's worker consumed it;
  an untracked scratch copy masked the loss for ten weeks, then vanished — recovery had to
  route through a surviving `.pyc`'s `co_filename` back into git history.)

### The mirror — landing a deferred capability stales its DEFERRAL CONTRACT

Retirement's mirror image, and just as blocking. When a change flips a case from *deferred*
to *implemented*, every docstring naming that case "raises / deferred / not yet supported"
becomes present-tense-FALSE, and the blast radius is the same three searches.

- The recurring half-cleanup signature: **the human-facing prose ledger gets rewritten and
  the machine-facing contracts do not** — the `@runtime_checkable` Protocol stub, the BASE
  class's default docstring that the next implementer inherits, the sibling class's own
  docstring, and public operators in files the diff never touched.
- Grep the deferred case's NAME and its prose forms across the package, then discriminate
  **by arm**: landing a matvec-transpose does not un-defer the transpose-SOLVE. Only the
  one flipped row changes; the still-genuine future seams stay.
- In a campaign-closing commit this is blocking — a `Closes #NN` trailer is internally
  inconsistent with a status the tree still tags deferred.

### Exception — keep a relinquished *fuller view* as a verification oracle

Retirement targets a *superseded* predecessor (same job, done worse). It does NOT target a
**fuller view of a concept that an optimization relinquished** (full field → rolling window;
full angular flux → moments; dense operator → factored form). The fuller view is a
**verification pathway** that pins the optimized path's reference — keep it.

- Make the keep-vs-retire decision EXPLICIT; never let a fuller view fall out silently or
  sit half-alive (orphaned-but-undeleted). It is either a wired, exercised oracle or it is
  retired.
- The oracle is NOT production-reachable and MUST share the optimized path's kernel (only
  the representation/storage differs), pinned by a permanent **end-to-end** equivalence test
  (`optimized ≡ fuller-view-oracle`), bit-identical or principled-equivalent per
  `vv-principles`.
- **Discriminator:** keep it only if it is a genuine structural reference. A "fuller view"
  that is the SAME math merely *procedurally rearranged*, AND already verified by a
  structurally-independent oracle elsewhere (MMS, closed form), is genuine redundancy —
  retire it. (`vv-principles` L11: procedural ≠ structural independence.) Worked: the 2-D
  rolling-window sweep kept `_sweep_2d_full_field` as a typed-`WavefrontFlux` oracle; the
  per-ℓ scattering kernel was retired as mere procedural rearrangement.
