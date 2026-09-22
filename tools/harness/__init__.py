"""The agent harness, generated from ``docs/development/``.

The docs are the source of the project's working knowledge — the rule cores,
the skill cores, the lessons index, each agent's definition, the on-boarding
page — and a harness (today Claude Code) reads a GENERATED, committed view of
them. The flow is one-way: the docs never depend on a harness's files (plan
``harness_context_budget.md``, K7, ruling 2026-09-19), so another harness
reads the same pages.

Modules, each answering one question:

* ``source``   — what is a source page? ``Kind``, ``Page``, discovery under
  ``docs/development/`` and the ``harness:`` front-matter block.
* ``links``    — how does a relative link move from the source's directory to
  the target's, and which anchors does MyST mint for a page?
* ``budget``   — what does a page cost, and what is a budget?
* ``render``   — the two shapes a target takes: a stamped whole file, or a
  block spliced between markers into a hand-maintained file.
* ``pipeline`` — generation for one harness; drift; orphans.
* ``targets``  — the ``Harness`` Protocol (``base``) and the implementations.
  The modules above name no harness and import none;
  ``tests/tools/test_harness_generator.py`` asserts it.
* ``__main__`` — ``python -m tools.harness [--check]``.
"""
