---
name: Stub-to-rich-narrative expansion sweeps for theory pages
description: Pattern for turning TODO-marker stub-grade theory pages into Cardinal-Rule-3 rich narratives. Three-page sweep (fn_method + singular_eigenfunction + sood_registry) doubled total LoC from 2252 to 4192 with zero Sphinx warnings; key disciplines documented below.
type: feedback
---

When a `docs/theory/<topic>.rst` page is marked Stub-grade with TODO
markers awaiting expansion, the maximum-effort task pattern is:

**Why:** Sphinx-as-the-LLM's-brain (Cardinal Rule 3) requires that
each page make a future reader an expert on the topic from the page
alone — full derivations, design rationale, what failed, gotchas,
literature with equation numbers, numerical evidence, cross-refs.
The TODO stubs that survived earlier session cycles are exactly the
content that's most expensive to recover from code, agent memory,
and primary literature each session — making them the highest-value
expansion target.

**How to apply (the playbook):**

1. **Read the entire stub page first** — every TODO has a brief
   inline pointer to its SymPy module + test gate. The brief is
   the prose seed; the SymPy docstring is the technical seed.
2. **Read the corresponding `origins/*_derivations.py` modules**
   (the algebra-of-record SymPy). Each `derive_*()` returns a dict
   with the LHS/RHS/diff/pass — the function docstring is usually
   80% of the prose you need.
3. **Read the closeout memos for any in-flight ERR-NNN entries**
   — these have the deepest context (mechanism, fingerprint,
   anti-pattern caught). For investigations that landed in the
   current session (e.g., ERR-036 Atkinson, ERR-037 μ=tanh, ERR-038
   Atalay floor) the memos at `.claude/agent-memory/numerics-investigator/`
   carry the rich-narrative fact density that the stub page should
   reflect.
4. **Replace each TODO note with its full prose** in this shape:
   - **SymPy derivation pointer** (`:func:`...).
   - **Test gate pointer** (`:func:`tests.derivations....`).
   - 1-2 paragraphs of math + design rationale + literature
     reference + the actual identity being verified.
   - For the load-bearing identities (e.g., the Atkinson F_k closed
     forms, the WM-72 q-formula, the Atalay parity flip), include
     the full derivation with intermediate steps.
5. **Add `:label:` blocks for every literature-transcribed equation**.
   Mark each with `.. vv-status: <label> documented` because they
   are definitions / transcriptions — not testable claims (the
   underlying identities are tested via SymPy ``simplify==0`` and
   the L1 numerical gates against published truth values, not
   against rendered LaTeX).
6. **The vv-status rationale comment** uses one of: `governing` (a
   transport equation defining the problem class — the entire
   solver IS the verification), `derivation` (an intermediate
   step in a verified chain), or for new categories, follow the
   peierls.rst pattern.
7. **Build Sphinx with `-W --keep-going` after every page** —
   warnings as errors. The capability matrix auto-regen +
   verification matrix auto-regen will register your new labels;
   the orphan count must stay 0 (every label is either
   `verifies("...")`-tested or `vv-status documented`).

**Common pitfalls during the sweep:**

- **Citation duplicate-definition warnings cross-document**. If a
  `[FooBar1979]_` citation is already defined in another theory
  page, redefining it in your expanded page produces a
  duplicate-citation warning. Use plain-text inline ("Foo & Bar
  1979") and add the reference as a regular bullet in the
  References section instead.
- **`:doc:` and `:ref:` paths**. Use relative `:doc:`fn_method``
  (no leading slash) when referencing sibling theory pages —
  absolute paths like `:doc:`/skills_of_record/algebra_of_record``
  fail because the destination doesn't exist as a Sphinx doc; the
  `algebra-of-record` skill is a `.claude/skills/` directory, not
  a Sphinx page.
- **Pygments lexer error on Unicode in code-blocks**. The `…`
  ellipsis (U+2026) is not a valid Python token; use `...` (three
  ASCII dots). Pygments emits "Lexing literal_block" warnings for
  any non-ASCII syntax it can't parse.
- **Math role with trailing space**. `:math:`c \leftrightarrow `
  Case's notation`` triggers a "phrase reference start-string
  without end-string" warning because the closing backtick after
  the space gets ambiguously parsed. Either close the math role
  cleanly (`:math:`c \leftrightarrow X``) or use a Unicode arrow
  in plain text ("c ↔ Case's notation").

**Quality scores from this sweep (1-5 scale):**

| Page | Derivation depth | Cross-refs | Numerical evidence | Failed approaches | Code traceability | Derivation source |
|------|:---:|:---:|:---:|:---:|:---:|:---:|
| fn_method.rst | 5 | 5 | 5 | 5 | 5 | 5 |
| singular_eigenfunction.rst | 5 | 5 | 5 | 5 | 5 | 5 |
| sood_registry.rst | 4 | 5 | 4 | 4 | 5 | n/a (registry, no derivation) |

The fn_method + singular_eigenfunction expansions hit 5/5 across
all dimensions because the closeout memos provided rich numerical
tables (ERR-036 convergence sweep, ERR-037 z_0 before/after by c,
ERR-038 1/d_crit scaling) and the SymPy modules provided the full
algebraic chain. sood_registry scored 4/4 on derivation/numerical
because it's a registry, not a solver — there's no "derivation"
per se to verify, only the schema invariants.

**Token-cost note:** the three pages totalled ~140K tokens of input
context (3 stub pages + 8 SymPy origin modules + 3 closeout memos +
ERR catalog entries) and ~70K tokens of generated output (4192
lines of expanded narrative). Budget ~3-4× the page's final word
count in input context, especially for V&V-rich pages where the
ERR-NNN cascade memos are the load-bearing source material.
