Development — how ORPHEUS is built, and the source of the agent harness
=========================================================================

This section is the **durable, harness-independent** home of the project's
working knowledge: the rules every contributor follows, the lessons and
founding cases behind them, the verification and elegance disciplines, the
agent workflows, and the git workflow. The ``.claude/`` directory that the
Claude Code harness reads is **generated from these pages** by
``tools/harness`` and never edited by hand; any other
harness (an API caller, a different agent runtime) reads the same pages.

Pages under ``rules/`` and ``skills/`` are the *cores* — what loads into an
agent's context at every session. Pages under ``evidence/`` carry the
founding cases, the surprise log and the lesson bodies, read on demand when
a core's link names them. See :doc:`harness` for what loads when and how to
add a rule, a skill or an agent.

New to the project, human or agent: :doc:`onboarding` first. It is the page
the harness's ``CLAUDE.md`` is generated from.

.. toctree::
   :maxdepth: 1
   :caption: Working with the repository

   onboarding
   git_workflow
   harness
   workflows
   lessons

.. toctree::
   :maxdepth: 1
   :caption: Rule cores (always-on)

   rules/cardinal
   rules/plan-authoring
   rules/coding-standards
   rules/instrument-doctrine
   rules/articulation
   rules/workflows
   rules/process-discipline
   rules/code-search
   rules/vv-testing

.. toctree::
   :maxdepth: 1
   :caption: Skill cores (loaded on demand or preloaded per agent)

   skills/vv-principles
   skills/coding-elegance
   skills/instrument-doctrine
   skills/retirement-audit

.. toctree::
   :maxdepth: 1
   :caption: Agent role blocks

   agents/index

.. toctree::
   :maxdepth: 1
   :caption: Evidence (founding cases, surprise log, lesson bodies)

   evidence/plan-authoring
   evidence/coding-standards
   evidence/process-discipline
   evidence/code-search
   evidence/lessons
   evidence/vv-anti-patterns
   evidence/test-design-modes
   evidence/coding-elegance
