Development Workflow
====================

ORPHEUS is developed in an **AI-heavy, solo-maintainer, Sphinx-as-brain**
style. This page formalizes the git workflow, commit conventions, and
branching discipline so that every session — human or AI — starts with
the same ground rules.

.. contents::
   :local:
   :depth: 2


Core principles
---------------

1. **``main`` is always green.** Always buildable, always Nexus-reloadable,
   tests passing. The Nexus MCP server reads the knowledge graph live —
   from wherever ``.nexus/config.toml`` declares it, which ``nexus config db``
   prints — so a broken ``main`` poisons every future agent session until
   it is fixed. Never commit directly to ``main``.

2. **History tells a story.** Commits within a branch should be
   granular enough to explain *why* each change happened. We preserve
   that story on merge by fast-forwarding, not squashing. Your future
   self (and every AI agent running ``git log``) reads the commit
   messages as documentation of *how the project got here*.

3. **No ceremony without value.** We do NOT use Git Flow
   (develop/release/hotfix branches), pull requests, code review, or
   protected branches. Those earn their keep in multi-developer teams.
   Solo AI-heavy development optimizes differently: the feedback loop
   is *human-agent*, not *human-human*. The one automation kept is the
   CI workflow (`Continuous integration`_ below), which runs the cheap
   gates on every push so that "``main`` is always green" is a claim
   someone reads.

4. **Sphinx is the contract.** Documentation changes and code changes
   are one commit, not two PRs. A feature is not DONE until the
   theory page, the test, and the Nexus graph all agree.


Branch model
------------

We use a **lightweight GitHub Flow** with a single permanent branch:

- ``main`` — the trunk; everything eventually lands here
- ``<type>/<topic>`` — short-lived feature branches, deleted on merge

Branch types follow the same vocabulary as commit prefixes (below):

.. csv-table::
   :header: Branch prefix, Used for, Example
   :widths: 15, 35, 50

   ``feature/``, new capability, ``feature/mms-solver``
   ``fix/``, bug or regression fix, ``fix/moc-attenuation-sign``
   ``docs/``, documentation-only work, ``docs/development-workflow``
   ``refactor/``, structural change with no behavior change, ``refactor/geometry-factories``
   ``test/``, test-only additions (rare; tests usually ride with feat/fix), ``test/err-catcher-retrofit``
   ``chore/``, tooling / dependencies / build config, ``chore/bump-nexus-0.9.0``

**Branch lifetime**: hours to days, not weeks. A branch that lives
longer than a week should be split or rebased against ``main``.

**Creating a branch**::

   git checkout main
   git pull --ff-only
   git checkout -b feature/<topic>


Commit message convention
-------------------------

We use **Conventional Commits** with optional scope. This makes
``git log --grep`` work as a queryable history, and unlocks automated
changelog generation if we ever want it.

**Format**::

   <type>(<scope>): <summary>

   <optional body — WHY, not WHAT>

   <optional trailer: Closes #NN, Co-Authored-By: ...>

**Types**:

.. csv-table::
   :header: Type, Meaning
   :widths: 15, 85

   ``feat``, "New user-visible capability (solver, method, data pipeline)"
   ``fix``, Bug fix for incorrect behavior
   ``docs``, Sphinx/RST/docstring/README changes only
   ``test``, Test-only additions or fixes
   ``refactor``, Structural change with no behavior change
   ``chore``, "Tooling, dependencies, build config, environment"
   ``perf``, Performance optimization without behavior change

**Scopes** (optional, one per commit): ``cp``, ``sn``, ``moc``, ``mc``,
``diffusion``, ``geometry``, ``data``, ``derivations``, ``kinetics``,
``fuel``, ``th``, ``numerics``, ``harness``, ``nexus``.

**Examples**::

   feat(cp): add interface current method for multi-cell lattice

   Closes #56.

   fix(moc): restore tau/sin_p factor in characteristic ODE

   ERR-007 reproduces under FM-07 (missing sine projection).
   L0 hand-calc test_attenuation_vacuum_source now catches it.

   docs(cp): document 9-case CP matrices in collision_probability.rst

   Closes #65.

   test(harness): retrofit @pytest.mark.catches for ERR-001..ERR-020

   refactor(geometry): extract zone-subdivision strategies

**Body** is optional for trivial changes. For non-trivial work,
explain **why** (what motivated the change, what alternatives were
considered, what failure mode this catches). The **what** is already
in the diff.

**Trailers**:

- ``Closes #NN`` auto-closes the GitHub issue on merge
- ``Co-Authored-By:`` lists AI co-authors when the work was
  agent-generated


Merging to ``main``
-------------------

Only one merge strategy: **fast-forward**.

::

   git checkout main
   git pull --ff-only
   git merge --ff-only <branch>
   git push origin main

If fast-forward fails (``main`` moved while you were working), rebase
the branch first::

   git checkout <branch>
   git rebase main
   # resolve conflicts, re-run tests
   git checkout main && git merge --ff-only <branch>

**Never merge-commit.** Never squash-merge. The linear history is
load-bearing: it lets ``git bisect`` pinpoint regressions to a single
commit, and it preserves the *chapter structure* that commits within a
branch naturally provide.

**After merging**, delete the branch (local + remote)::

   git branch -d <branch>
   git push origin --delete <branch>


Continuous integration
----------------------

One GitHub Actions workflow, ``.github/workflows/gates.yml``, runs on every
push and every pull request (landed 2026-09-21). It is the cheap,
deterministic half of the project's acceptance, in this order:

1. ``python -m tools.harness --check`` — the harness view (``CLAUDE.md``,
   ``.claude/``) is generated from ``docs/development/`` and current. A hand
   edit inside a generated file, or a source edited without regeneration,
   reds here.
2. ``python -O -m pytest tests/tools/test_harness_generator.py
   tests/test_harness_generated.py tests/test_layer_imports.py`` — the
   generator's own tests, the generated-tree gate and the layer contract
   (417 tests; [M] 2026-09-21: 6.4 s locally).
3. ``python -O -m pytest tests/test_pyright_ratchet.py`` — the type gate as
   the #226 ratchet: pyright 1.1.410 (installed globally in the job) over
   ``orpheus/``, its per-module error counts compared with the committed
   baseline ``tests/_harness/pyright_baseline.json`` (0 errors). Pyright
   reads ``pyproject.toml``'s ``[tool.pyright]``, which names ``.venv``; the
   workflow installs into a ``.venv`` of its own for that reason, so CI and a
   developer's machine read one configuration. The scope is ``orpheus/``
   only: ``tests/`` and ``tools/`` are outside the gate today (about 1 400
   pre-existing errors in ``tests/``, #226's burn-down), and a bare
   ``pyright`` over the whole tree is not a gate anywhere (it walks
   ``scratch/`` and every worktree; [M] 2026-09-06 it exhausts node's heap
   locally, and [M] 2026-09-21 the first CI run read 1 674 errors from it).
4. ``sphinx-build -E -W --keep-going docs docs/_build/html`` — the strict
   build, which also regenerates every generated file and the Nexus graph
   ([M] 2026-09-21: 77 s locally).
5. ``git diff --exit-code`` — the build changed no tracked file, so every
   generated artefact committed to the tree was current.

The full pytest suite is NOT in CI. It is serial and takes over ninety
minutes per tier, and it stays the local gate run before a merge (the
canonical invocation is ``python -O -m pytest``). The workflow exists so
that the ``process-discipline`` clause "after pushing, look at CI" has an
instrument to read. After ``git push origin main``::

   gh run list --branch main --limit 3
   gh run watch            # follows the latest run to completion

A red run on ``main`` is fixed before anything else lands, and a run that
was already red before your push is baselined first (the rule's baseline
clause). The workflow landed with its first red: a deliberate hand edit
inside a generated rule, pushed on the workflow's branch and reverted, so
that the green run after it is a baseline and not an untested instrument
(instrument doctrine X1).


Release tagging
---------------

Tag ``main`` at meaningful V&V milestones::

   git tag -a v0.2 -m "v0.2: full SN + CP + MOC verification ladder"
   git push origin v0.2

Tags let Nexus snapshot "the state of the graph when equation X was
first verified" and give humans a referenceable point for teaching
material. No strict semver yet — ORPHEUS is pre-1.0, so minor-version
bumps are meaningful-milestone markers, not API-stability promises.


When to branch and when not to
------------------------------

**Branch** for anything that:

- touches more than one file, OR
- takes more than one commit, OR
- might not work on the first try and needs rollback safety, OR
- spans more than one conversational turn

**Don't branch** for:

- typo fixes
- adding a single test
- updating a single reference in an existing doc
- emergency ``main`` fixes that are smaller than a branch roundtrip

These can go straight on ``main`` — but they are rare. When in
doubt, branch. The cost of a branch is one ``git checkout``; the
cost of a polluted ``main`` is a broken Nexus graph for every
subsequent session.


AI-agent workflow
-----------------

When a Claude Code session starts a non-trivial task:

1. **Read** the ``CLAUDE.md`` session-start protocol at the repo root.
2. **Check issues** relevant to the module.
3. **Plan** the work (use plan mode for 3+ steps).
4. **Branch** before the first edit: ``git checkout -b <type>/<topic>``.
5. **Implement** — code + tests + docs + Sphinx, one commit per
   logical step.
6. **Verify** — run ``python -O -m pytest``, rebuild Sphinx strictly
   (``sphinx-build -E -W --keep-going docs docs/_build/html``), run
   ``python -m tests._harness.audit``.
7. **Merge** — ``git merge --ff-only`` to ``main``, push, delete branch,
   then read the CI run (``gh run watch``; `Continuous integration`_).
8. **Close** linked issues via commit trailers or explicit
   ``gh issue close``.

Sub-agents inherit this workflow implicitly: the CLAUDE.md session
protocol anchors them to the same rules, and the Nexus graph indexes
this page so ``nexus-exploring`` queries surface it automatically.


Quick command reference
-----------------------

**Start a branch**::

   git checkout main && git pull --ff-only && git checkout -b feature/<topic>

**Commit with prefix** — through ``-F`` from a quoted heredoc, never ``-m``: zsh
command-substitutes backticks inside ``-m "…"`` and the words vanish silently; read
the message back with ``git log -1 --format=%B``::

   git commit -F - <<'MSG'
   feat(cp): add interface current method
   MSG

**Merge and clean up**::

   git checkout main && git merge --ff-only feature/<topic> && git push origin main
   git branch -d feature/<topic>
   git push origin --delete feature/<topic>

**Query history**::

   git log --grep '^feat'        # every feature
   git log --grep '^fix'         # every bug fix
   git log --grep '^docs'        # every doc-only change
   git log --grep '(cp)'         # every CP-module change

**Verify before merge**::

   python -O -m pytest tests/ -q
   sphinx-build -E -W --keep-going docs docs/_build/html
   python -m tests._harness.audit

**Read the CI after a push**::

   gh run list --branch main --limit 3
   gh run watch


Related pages
-------------

- :doc:`/theory/verification/harness` — V&V test harness conventions
- :doc:`/theory/verification/index` — the verification part, with the
  auto-generated V&V matrix
- ``CLAUDE.md`` — session-start protocol and cardinal rules
- ``docs/theory/verification/error_catalog.rst`` — caught-bug publication artifact
