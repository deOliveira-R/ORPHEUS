---
name: Auto-generated Sphinx tables — registry as single source of truth
description: Pattern for replacing hand-written tables with registry-driven includes (generator tool + capability_rows metadata function + conf.py builder-inited hook)
type: feedback
---

When a hand-written Sphinx table enumerates registry data (e.g. the
Peierls capability matrix that mirrors ``continuous_all()`` filtered
to an operator form), ship three pieces rather than keep the table
hand-edited: (1) a metadata-only function on the registry side
(e.g. ``peierls_cases.capability_rows()``) that returns the row
dicts WITHOUT running any expensive solves; (2) a generator tool in
``tools/verification/generate_*_matrix.py`` that renders the rows as
an RST ``.. list-table::`` include file under ``docs/theory/`` (or
``docs/_generated/`` when the theory page has no sibling includes);
(3) a ``setup(app) -> app.connect("builder-inited", ...)`` hook in
``docs/conf.py`` that invokes the tool before Sphinx collects
sources. The existing ``tools/verification/generate_matrix.py`` is
the pattern template — its conf.py hook is already wired so you
just append a second connect().

**Why:** The peierls capability matrix was last drifting when the
Issue #104 commit added 2G cyl/sph references and the hand-written
table lagged by 3 commits; landing an auto-generator eliminates
drift (consolidation task T2.1). The trap is that calling
``continuous_cases()`` to get the metadata costs O(minutes) because
each ``ContinuousReferenceSolution`` construction runs a full
eigenvalue solve — that is unacceptable per-build cost. A sibling
metadata-only function (``capability_rows()``) that lists names,
shapes, and closure labels by inspecting the static loop structure
and imported tolerance tables is the cheap substitute; a test-suite
cross-check that the metadata list agrees with ``continuous_cases()``
row-for-row keeps them from drifting apart.

**How to apply:** Before hand-updating a registry-backed Sphinx
table, check for a generator under ``tools/verification/``. If the
table has 5+ rows or any parameter sweep, it is a candidate for
auto-generation. Always use ``.. list-table::`` (never csv-table)
for the generated output because many cells contain LaTeX with
commas. Keep Class A / Class B closure-per-shape split tables
hand-written — those document structural claims per shape rather
than shipped rows per reference, so they drift slowly.

**RST pitfall re-observed 2026-04-23:** A ``.. _label:`` anchor
must precede a **section title** (or a named directive like
``.. figure::``) to resolve as ``:ref:`label```. Attaching it to a
bold-text paragraph silently creates the label but every
cross-reference fails with "Failed to create a cross reference.
A title or caption not found". Fix: wrap the content in a titled
sub-section. Feedback_peierls_docs.md's title-underline code-point
rule still applies.

**Sphinx section marker depth in peierls_nystrom.rst (2026-04-23):**
``====`` = h1, ``----`` = h2, ``~~~~`` = h3, ``^^^^`` = h4. Skipping
a level (e.g. inserting ``^^^^`` directly under ``----``) raises
ERROR: Inconsistent title style: skip from level 3 to 5. Always match
the file's established ladder — grep the file for a prior usage of
``~`` or ``^`` underlines and copy the pattern.

**Duplicate labels across files:** If a label name already exists
elsewhere (e.g. an orphan ``.. _peierls-scattering-convention:`` in
collision_probability.rst that was never attached to a title), and
you add the same label to a new file, Sphinx emits "duplicate
label" at build time. Consolidation strategy: keep the label on the
page that is meant to be the canonical source; convert all other
occurrences to ``:ref:`label``` cross-references (not label
definitions).
