"""Hand-written mutation batteries — the instruments that prove gate teeth.

A *battery* here is a small set of **hand-authored, named, hazard-specific**
mutations, each the plausible transcription of a real bug, run as an in-process
pytest plugin. See ``README.md`` for what this is, what it is NOT
(``tests/_mutation/`` is a different mechanism), and how to run one.

Nothing in this package is collected by pytest: no module matches ``test_*``.
The modules are *plugins*, loaded explicitly with ``-p``.
"""
