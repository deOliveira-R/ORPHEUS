.. _theory-references-bibliography:

============
Bibliography
============

This page is the ORPHEUS theory corpus's **single citation home**. Every
``:cite:`` role across ``docs/theory/`` resolves to an entry on this page
and links here.

The bibliographic records themselves live in ``docs/refs.bib`` — the single
source of truth, maintained upstream in Zotero (manual Better-BibTeX export).
Citation labels render as the BibTeX keys verbatim (``bibtex_default_style``
is set to ``keylabel`` in :file:`docs/conf.py`), so a citation of the key
``Hebert2009`` displays as ``[Hebert2009]`` — character-identical to the
docutils bracket labels this page replaced.

The former per-page "References" sections (docutils ``.. [Key]`` definition
blocks scattered across the theory pages) were retired in favour of this
centralised bibliography during the ``sphinxcontrib-bibtex`` migration
(Issue #231, Phase G2). To add, amend, or correct a reference, edit
``docs/refs.bib`` (or its Zotero source) — **never** a theory page. Only
works actually cited somewhere in the corpus appear below.

.. bibliography::
