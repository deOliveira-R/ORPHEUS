r"""A citation of a published source: its bibliography key and where in it.

A :class:`Citation` names a work by its key in the project's bibliography,
``docs/refs.bib`` (the file ``sphinxcontrib.bibtex`` renders, and the same
key a docstring cites with ``:cite:``), and optionally a place inside the
work: a problem, a table, an equation or a page. It is the typed form of the
free-text citations the reference registries used to carry.

Each citation goes where its claim lives (the reference-solution campaign,
``.claude/plans/reference_cache.md``). The source that defined a benchmark
problem is cited by the problem; the source of a published value is cited by
that value. The literature of a method (the equation a generator implements)
is cited in the generator's docstring and is not copied into data.

The key is parsed at construction against the grammar every key in
``docs/refs.bib`` follows, a letter then letters, digits or underscores
(79 of 79 keys, 2026-09-29). Whether the key exists in the file is a gate
(``tests/gates/derivations/test_registry_citations_resolve.py``), not a
runtime check: ``docs/`` does not ship with the package.

The edition of a work is part of its key. Sood, Forster and Parsons' 1999
report and their 2003 journal paper are two entries, ``SoodLA13511_1999``
and ``SoodForsterParsons2003``, whose tables, equations and references are
numbered differently, so a locator is meaningful only with its key.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

_BIBKEY = re.compile(r"[A-Za-z][A-Za-z0-9_]*")


@dataclass(frozen=True)
class Citation:
    """A published work, by its ``docs/refs.bib`` key, and a place in it.

    Parameters
    ----------
    bibkey : str
        The work's key in ``docs/refs.bib``: a letter, then letters, digits
        or underscores.
    locator : str or None
        Where in the work, in the work's own numbering (``"problem 4"``,
        ``"Table 10"``, ``"p. 71"``); ``None`` cites the work as a whole.

    Raises
    ------
    ValueError
        If the key does not follow the grammar, or the locator is blank.
    """

    bibkey: str
    locator: str | None = None

    def __post_init__(self) -> None:
        if not _BIBKEY.fullmatch(self.bibkey):
            raise ValueError(
                f"Citation bibkey {self.bibkey!r} is not a refs.bib key: a key is a "
                "letter followed by letters, digits or underscores"
            )
        if self.locator is not None and not self.locator.strip():
            raise ValueError(
                f"Citation locator for {self.bibkey!r} is blank: give a place in "
                "the work, or None to cite the work as a whole"
            )


__all__ = ["Citation"]
