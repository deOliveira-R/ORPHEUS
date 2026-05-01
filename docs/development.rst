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
   tests passing. The Nexus MCP server reads from ``docs/_build/html/_nexus/graph.db``
   live, so a broken ``main`` poisons every future agent session until
   it is fixed. Never commit directly to ``main``.

2. **History tells a story.** Commits within a branch should be
   granular enough to explain *why* each change happened. We preserve
   that story on merge by fast-forwarding, not squashing. Your future
   self (and every AI agent running ``git log``) reads the commit
   messages as documentation of *how the project got here*.

3. **No ceremony without value.** We do NOT use Git Flow
   (develop/release/hotfix branches), pull requests, code review, or
   protected branches. Those earn their keep in multi-developer teams
   with CI gates. Solo AI-heavy development optimizes differently:
   the feedback loop is *human-agent*, not *human-human*.

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
6. **Verify** — run ``pytest``, rebuild Sphinx, run
   ``python -m tests._harness.audit``.
7. **Merge** — ``git merge --ff-only`` to ``main``, push, delete branch.
8. **Close** linked issues via commit trailers or explicit
   ``gh issue close``.

Sub-agents inherit this workflow implicitly: the CLAUDE.md session
protocol anchors them to the same rules, and the Nexus graph indexes
this page so ``nexus-exploring`` queries surface it automatically.


Quick command reference
-----------------------

**Start a branch**::

   git checkout main && git pull --ff-only && git checkout -b feature/<topic>

**Commit with prefix**::

   git commit -m "feat(cp): add interface current method"

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

   pytest tests/ -q
   sphinx-build -q docs docs/_build/html
   python -m tests._harness.audit


Related pages
-------------

- :doc:`testing/architecture` — V&V test harness conventions
- :doc:`verification/index` — auto-generated V&V matrix
- ``CLAUDE.md`` — session-start protocol and cardinal rules
- ``.claude/skills/vv-principles/error_catalog.md`` — caught-bug publication artifact
