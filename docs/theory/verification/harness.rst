.. vv-audit: skip-file
   (this page teaches the ``:label:`` / ``.. vv-status:`` syntax with
   verbatim examples; the marker keeps the audit's theory scan from
   reading those examples as data — see "The audit CLI" below)

Test-Harness Architecture
=========================

.. contents::
   :local:
   :depth: 2

Motivation
----------

ORPHEUS operates a four-level physics-verification ladder (L0..L3)
plus an orthogonal ``foundation`` bucket for software-invariant
tests that don't correspond to a physics equation. The ladder
itself — what each rung proves, what it deliberately does not
prove, and the necessity chain between rungs — is defined
normatively at :ref:`vv-level-ladder` in :doc:`principles`;
``foundation`` is not on the ladder — see
:ref:`vv-foundation-tests` for the taxonomy. This page owns the
*harness* side of the contract: how a test declares its rung, how
the declarations are audited, and how the declared facts reach the
verification matrix and the Nexus knowledge graph.

One ladder rule is worth restating wherever tests are authored:
1-group eigenvalue tests are **degenerate** for transport
verification — :math:`k = \nu\Sigma_f / \Sigma_a` regardless of
flux shape — so always demand ≥2 groups for any test claiming to
verify transport. The canonical statement and the ERR-006 war
story live at :ref:`verification-1g-degeneracy`.

Design principles
-----------------

1. **One source of truth per fact.** The V&V level a test belongs to,
   the equation labels it verifies, and the failure modes it catches
   are declared *once*, on the test itself, as pytest markers. The
   audit tool, the Sphinx verification-matrix page, and the Nexus
   knowledge graph all consume the same declaration.

2. **No new DSL.** The convention is vanilla ``pytest.mark.*`` —
   one spelling of every declaration, tree-wide.

3. **Inherit from the reference case when possible.** Tests that pull
   analytical values via :func:`ref` inherit the V&V level from the
   underlying :class:`~orpheus.derivations.common.verification_case.VerificationCase`.
   Case metadata is populated once (in
   :mod:`orpheus.derivations`) and every consuming test is tagged
   automatically by the conftest hook.

4. **Nexus-native traceability.** Tests reference equations via
   ``:math:`label``` docstring roles. sphinxcontrib-nexus ≥ 0.6.0
   converts those into graph edges from the test node to the
   corresponding ``math:equation:*`` node, so
   ``verification_coverage`` answers "which tests verify X" directly
   from the declared markup.

5. **Central audit, not per-file checks.** A single command —
   ``python -m tests._harness.audit`` — produces the full V&V matrix,
   lists orphan equations, and cross-checks
   ``.claude/skills/vv-principles/error_catalog.md``. No scattered assertions.

6. **Enforcement mode.** Every test in ``tests/`` carries a level
   tag — physics (``l0``..``l3``) or ``foundation``. The audit tool
   surfaces every untagged test and ``--strict`` exits non-zero on
   any gap, so new tests cannot slip in untagged. The "unmarked
   accumulates in its own bucket" stance from the initial migration
   is behind us; foundation finally gave the non-physics tests a
   home (see :ref:`vv-foundation-tests`).

7. **Type-error ratchet** (issue #226). The package carries a large
   pre-existing pyright error surface; until the per-module burn-down
   reaches zero, the enforceable invariant is monotonicity.
   ``tests/test_pyright_ratchet.py`` (``foundation`` + ``slow``,
   skips without a host pyright) compares live per-module error
   counts against ``tests/_harness/pyright_baseline.json`` and fails
   in BOTH directions — an increase is a regression, a decrease must
   be locked in via
   ``python -m tests._harness.pyright_ratchet --update`` so
   improvements can't silently erode. The single source for the
   counting is :mod:`tests._harness.pyright_ratchet`; the baseline
   records the pyright version because counts move across pyright
   releases without code changes.

Authoring a test
----------------

Raw ``pytest.mark.*`` is the ONE convention in the ORPHEUS codebase
(every test file uses it). If you are writing a new test, raw
markers are the path — there is no alternative spelling.

Raw ``pytest.mark.*`` decorators
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the convention used by every test file in the tree.

.. code-block:: python

   import pytest


   @pytest.mark.l0
   @pytest.mark.verifies("transport-cartesian")
   @pytest.mark.catches("FM-07", "ERR-003")
   class TestSingleTrackAttenuation:
       """L0: Verify :math:`transport-cartesian` for a characteristic track.

       For a pure-absorber slab of thickness L with vacuum inlet, the
       exit flux is :math:`\\psi_\\text{out} = \\psi_\\text{in} \\cdot
       e^{-\\Sigma_t L / \\mu}`. Catches sign flips in the
       :math:`\\Delta\\psi` update (FM-07) and missing :math:`\\tau`
       normalization.
       """

       def test_attenuation_vacuum_source(self): ...
       def test_attenuation_equilibrium(self): ...

The class-level ``@pytest.mark.l0`` and ``@pytest.mark.verifies(...)``
cascade to every ``test_*`` method inside the class. The
docstring's ``:math:`transport-cartesian``` role is picked up by
sphinxcontrib-nexus and written as a graph edge from the test node
to ``math:equation:transport-cartesian`` on the next ``sphinx-build``.

For file-level application (the most common shape in the repo —
see ``test_cp_verification.py`` or ``test_homogeneous.py``), use
``pytestmark`` at module scope:

.. code-block:: python

   pytestmark = [pytest.mark.l1, pytest.mark.verifies(
       "collision-rate", "p-inf", "matrix-eigenvalue", "mg-balance",
   )]

Foundation tests use ``@pytest.mark.foundation`` instead of an
``lN`` marker and never declare ``verifies(...)``:

.. code-block:: python

   pytestmark = pytest.mark.foundation  # file-level, test_geometry.py

(A ``verify.lN(...)``/``vv_cases(...)`` sugar layer existed until
2026-07 but was retired with zero consumers — a second spelling of
the same declaration is a twin path, not an ergonomic win.) Tests
that parametrize directly over case objects
(``@pytest.mark.parametrize("case", [...])`` where each object has a
``vv_level``) inherit their level through the conftest hook exactly
like the ``ref()`` shape below.

Inheriting through ``ref()``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Tests that use the ``ref`` fixture inherit the V&V level and
equation labels from the underlying case automatically. **No explicit
marker needed** when a test is a thin wrapper around a single case:

.. code-block:: python

   def test_homogeneous_1g_eigenvalue(ref):
       case = ref("homo_1eg")   # vv_level="L1", equation_labels=("matrix-eigenvalue",)
       ...

At collection time the conftest hook reads ``case.vv_level`` and
``case.equation_labels``, stamps the test with ``@pytest.mark.l1``
and ``@pytest.mark.verifies("matrix-eigenvalue")``, and records the
provenance as ``level_source="case"`` in the registry.

Precedence order (most specific wins)
-------------------------------------

The conftest hook applies V&V level markers in this order:

1. Explicit marker on the test itself (``@pytest.mark.lN`` or
   file-level ``pytestmark``).
2. Class name matching ``TestL<N>Foo`` (legacy convention, still
   honored).
3. Function name matching ``test_l<N>_*``.
4. Inherited from :class:`VerificationCase` via a parametrized
   ``case`` argument.
5. **Unmarked** — recorded in the registry with ``level=None`` so the
   audit tool can surface it.

Conflicts between different ``lN`` markers on the same test are
resolved deterministically (highest level wins) and a warning is
emitted so duplicates surface immediately. The ``foundation``
marker sorts below every ``l<N>`` in the tiebreak, so a test
accidentally carrying both ``l1`` and ``foundation`` resolves to
``L1`` (the stronger physics claim) and the foundation marker is
surfaced as a conflict. See :ref:`vv-foundation-tests` for why
foundation is orthogonal to the physics ladder and why physics
always wins the tiebreak.

The audit CLI
-------------

.. code-block:: console

   $ python -m tests._harness.audit
   ========================================================================
   ORPHEUS V&V Test Audit
   ========================================================================
   Total tests collected: NNNN

   By V&V level:
     L0           ...
     L1           ...
     L2           ...
     L3           ...
     foundation   ...
     unmarked     ...

   By tagging source:
     explicit / class-name / func-name / case / unmarked

   Module × level grid:
     module                                 L0   L1   L2   L3   FD   ??
     ------------------------------------------------------------------
     cp/test_properties                     ...
     sn/test_properties                     ...

   Equation coverage:
     <label>                                   N test(s)

   Orphan equations (N of M testable theory labels have zero test
   coverage; K labels are .. vv-status: documented and excluded from
   the orphan gate):
     ...

   error_catalog.md ERR coverage (N/N entries have a catching test):

The output shape above is illustrative — every count drifts with each
commit, so run the CLI for current numbers (the auto-generated
verification matrix page carries the same data, refreshed on every
Sphinx build). The ``FD`` column counts
:ref:`foundation tests <vv-foundation-tests>`. Documented-only labels
are marked via the ``.. vv-status: <label> documented`` sentinel
described in :ref:`vv-status-documented`.

The tool runs ``pytest --collect-only`` under the hood so the
:data:`tests._harness.registry.TEST_REGISTRY` is populated, then
queries it. No test code is executed.

Flags:

``--json``
    Machine-readable output (full registry dump plus orphan /
    documented / phantom / ERR-coverage sets).
``--untagged``
    List only tests with ``level=None``. Should return an empty list
    under normal operation; non-empty output means new tests were
    added without a V&V tag (``l0``..``l3`` or ``foundation``).
``--gaps``
    List orphan equations (labels under ``docs/theory/**/*.rst`` with
    zero verifying tests, excluding documented-sentinel labels),
    phantom verifies-targets (tests naming a ``:label:`` that exists
    nowhere under ``docs/``), and ``ERR-NNN`` entries in
    ``.claude/skills/vv-principles/error_catalog.md`` with no
    catching test.
``--strict``
    Exit 1 if **any** of three gates trip:

    1. untagged tests exist (no ``l0``..``l3`` / ``foundation`` marker),
    2. orphan equations exist (theory labels with no ``verifies(...)``
       decorator pointing at them, ignoring documented-sentinel labels),
    3. phantom verifies-targets exist (a label rename/removal not
       migrated into its tests).

    Missing ERR catchers are reported but not strict-gated. The
    orphan backlog is being classified per-label under the V&V-part
    consolidation (task #10); until it lands, ``--strict`` is
    informational and the harness is run by hand before every merge
    (there is no CI).

    The ``--strict`` gate ignores any theory label that is marked
    :ref:`vv-status-documented` — those are deliberately excluded
    from the orphan set because they cannot or should not be paired
    with a test. A real gap (implemented-but-untested equation)
    still fires the gate.

    Independent of ``--strict``, **sentinel violations are a hard
    error on every invocation** (exit 2, before collection): an
    unknown vv-status word, a sentinel whose label is missing from
    its own file, or a malformed sentinel line each abort the audit —
    and therefore fail the Sphinx build that regenerates the matrix
    (fatal under ``-W``).

Scan-exempt files — the ``vv-audit: skip-file`` marker
------------------------------------------------------

A file under the theory tree can opt out of the label/sentinel scan
with a column-0 comment anywhere in its source::

   .. vv-audit: skip-file

The scanner is line-based — it cannot tell a literal-block *teaching
example* of the ``:label:`` / ``.. vv-status:`` syntax from the real
thing — so exactly two pages carry the marker: **this page** (its
sentinel and label blocks are verbatim syntax examples, not
declarations) and the **generated matrix page** (the generator emits
the marker; the page's label mentions are prose about the census,
not members of it). The audit reports every skipped file in all of
its output modes and the matrix lists them in its "Scan-exempt
files" section — the exclusion is always visible, never silent.

Never mark a real theory page: hiding genuine equations from the
orphan gate is exactly the silent-drop failure the fail-loud
sentinel schema exists to prevent. If a real page's example code
ever trips the scanner, the example is the thing to restructure.

.. _vv-foundation-tests:

Foundation tests — software invariants outside the L0..L3 ladder
-----------------------------------------------------------------

The L0..L3 ladder (:ref:`vv-level-ladder`) is organized around
Cardinal Rule 4 — "Are we solving the equations right?" Each rung
assumes there is a **physics equation** in a Sphinx theory page
being verified: L0 checks a
hand-calculation of a single term, L1 asserts measured convergence
order against an analytical or manufactured solution, L2 proves
multi-group heterogeneous consistency, and L3 compares against
experiment. A test that doesn't verify a labelled theory equation
has no natural place on this ladder.

But some tests exist that are **not** about physics:

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Example
     - What it verifies
   * - ``test_geometry::test_cartesian_single_cell``
     - ``compute_volumes_1d`` returns the right cell volume for
       known edges — a data-structure contract of the
       :class:`~orpheus.geometry.Mesh1D` factory.
   * - ``test_structured_geometry::test_equal_volume_{cylindrical,spherical}_invariant``
     - Every cell in an equal-volume zone has bit-identical volume
       by construction (the algebraic invariant that caught
       ERR-020). Not a physics claim — a round-trip
       correctness property of ``Mesh1D.from_geometry`` equal-volume
       subdivision.
   * - ``test_geometry::TestPWRPin2D``
     - ``StructuredGeometry.wigner_seitz_pin_cell`` returns the
       correct cell radius from the Wigner-Seitz identity. Geometric
       primitive, not transport physics.
   * - ``test_geometry::TestMesh1D``
     - ``Mesh1D`` instances are frozen and reject invalid inputs.
       Language-level contract, not a reactor-physics equation.

These tests **must** exist — every downstream solver depends on
them — but there is no theory label for "volumes are computed
correctly" or "the factory rejects non-monotone edges." Calling
them L0 is a category error: L0 is term verification of a physics
equation, and these tests don't have a physics equation to verify.

The ``foundation`` marker exists for exactly this case.

**When to use** ``@pytest.mark.foundation``:

- The test verifies a **software invariant**: data-structure
  contract, numerical primitive, factory output, immutability
  guard, input validation, algebraic identity of a pre-physics
  building block.
- There is **no** ``:label:`` in any ``docs/theory/**/*.rst`` page
  that corresponds to what the test is checking. Foundation tests
  **never** carry ``@pytest.mark.verifies(...)`` — if they did,
  they would belong on the L0..L3 ladder instead.
- A failure means the code is broken **as software**, not that
  the physics is wrong.

**When NOT to use** it:

- The test does have a natural theory-page label. Use L0..L3
  instead — the physics ladder is the stronger claim.
- You aren't sure what level a test should be. The anti-pattern
  is "I don't know what level this is, so I'll call it
  foundation." If in doubt, read the test's docstring: if it
  reads like "term X of equation Y matches a hand calculation,"
  it's L0; if it reads like "this data structure satisfies
  invariant Z," it's foundation.
- The test is testing a reference implementation or derivation
  helper. Those belong in the derivation scripts' own tests, not
  the main suite.

**Interaction with other markers.** ``foundation`` is orthogonal
to the physics ladder. If two conflicting markers are applied
(e.g. ``l1`` and ``foundation`` on the same test), the physics
level wins — see the tiebreak rule in ``_existing_level`` of
``tests/conftest.py``. The ``catches("ERR-NNN")`` decorator is
orthogonal to the level bucket: a foundation test can absolutely
be the catcher for an ERR entry (ERR-020 is caught by
``TestZoneSubdivision``, which is foundation).

**Audit reporting.** ``python -m tests._harness.audit`` reports
foundation tests in their own row of the level breakdown and
their own ``FD`` column of the module × level grid, separate
from L0..L3. They do not affect the orphan-equation gate
(foundation tests carry no ``verifies(...)``, so they cannot
close an orphan) but they do satisfy the ``--strict`` mode's
untagged-tests gate — a foundation marker is a valid tag.

**Selection at runtime.**

.. code-block:: console

   pytest -m foundation              # only foundation tests
   pytest -m "not foundation"        # only physics V&V
   pytest -m "l0 or foundation"      # L0 + foundation (fast)


.. _vv-status-documented:

Documented-only equations (``:vv-status:``)
-------------------------------------------

Not every equation in ``docs/theory/**/*.rst`` can or should carry a
verifying test. Three cases come up in practice:

1. **Pure definitional labels.** ``boltzmann``, ``transport-equation``,
   ``balance-general`` — these name the governing equation or a
   mathematical identity. They have no single "implementing
   function" to test against; the entire transport solver *is* the
   verification, and the individual labelled test is exercised by
   downstream equations like ``matrix-eigenvalue`` and ``mg-balance``.
2. **Not-yet-implemented modules.** A theory page may document the
   full equation set of a module whose Python port does not yet
   exist (the TH / fuel-behaviour / reactor-kinetics modules are
   currently in this state — they live in the ``docs/theory``
   narrative but not in the ``orpheus/`` package tree). A
   documented-but-not-implemented equation is a work-in-progress
   marker, not a V&V gap.
3. **Deliberately deferred tests.** When writing the catching test
   requires infrastructure that does not exist yet (a new analytical
   reference, a missing fixture, a dependency to land in a separate
   issue), marking the label as documented-only is the escape hatch.
   This should be rare and each case should reference a tracking
   issue in the RST comment.

The V&V harness recognises these via a plain RST comment of the
form

.. code-block:: rst

   .. math::
      :label: boltzmann

      \partial_t \psi + \Omega \cdot \nabla \psi = S - \Sigma_t \psi

   .. vv-status: boltzmann documented

Because the line starts with ``.. `` followed by text that is **not**
a registered Sphinx directive, Sphinx silently strips it from the
rendered output — the sentinel lives only in the source file. The
audit CLI parses these comments and excludes the named labels from
the ``Orphan equations`` count and the ``--strict`` gate.

Rules (fail-loud since the 2026-07 single-status ruling, task #10):

- The ``vv-status:`` comment must appear in the **same RST file** as
  the ``:label:`` it refers to — a sentinel is point-of-use metadata.
  The audit **enforces** this: a sentinel whose label lives in a
  different file is a violation (the message says which file to move
  it to).
- ``documented`` is the **only** status, by design. ``tested`` /
  ``verified`` are *derived* facts — the matrix computes them from
  ``@pytest.mark.verifies`` declarations — so a hand-written coverage
  claim would be a second source of truth that can silently lie. Any
  other status word is a hard audit error, not a no-op.
- A sentinel naming a label that exists nowhere (a typo, or a label
  renamed without migrating the sentinel) is a hard audit error.
- Every violation aborts the audit with exit 2 **before** collection,
  which fails the matrix regeneration and therefore the Sphinx build
  (fatal under ``-W``) — invalid V&V metadata can never sit silently
  in the tree.
- Do not use the documented sentinel to paper over a genuine
  gap. "The test is hard to write" is not a justification;
  "the code does not exist yet" or "this is a definitional label"
  are. If in doubt, open an issue referencing the label.


Selecting tests at runtime
--------------------------

The standard pytest marker expressions apply:

.. code-block:: console

   pytest -m l0                       # only L0 term verification
   pytest -m "l1 and not slow"        # fast L1 checks
   pytest -m "l2 or l3"               # integration + validation
   pytest -m foundation               # only foundation tests (software invariants)
   pytest -m "not foundation"         # only physics V&V
   pytest -m "l0 or foundation"       # L0 + foundation (fast; excludes eigenvalue runs)
   pytest -m "verifies and not slow"  # any test with an equation label

Since ``verifies`` and ``catches`` are pytest marks with arguments,
``pytest -m "verifies"`` selects every test carrying any such mark.
Filtering by a specific label requires the audit tool (pytest's mark
expression language doesn't parse marker arguments).

``tests/_harness`` package layout
---------------------------------

.. code-block:: text

   tests/_harness/
       __init__.py            # re-exports TEST_REGISTRY, TestMetadata
       registry.py            # TestMetadata dataclass + TEST_REGISTRY dict
       audit.py               # python -m tests._harness.audit
       predicates.py          # two-axis inverse/adjoint operator contract
       pyright_ratchet.py     # #226 pyright error-count monotonicity gate
       pyright_baseline.json  # the ratchet's committed baseline
       xs.py                  # shared cross-section builders (re-exports)
       meshes.py              # (stub) shared mesh/geometry builders

``xs.py`` re-exports the canonical cross-section helpers from
``orpheus.derivations.common.xs_library`` (``make_mixture``, ``get_mixture``,
``get_xs``, ``get_materials``, ``validate_all``) so tests can import
them from a single stable path. ``meshes.py`` is currently an empty
placeholder — the shared ``_ws_mesh``, ``_homogeneous_ws_mesh``, and
related helpers are still duplicated across ``test_moc_verification.py``,
``test_cp_verification.py``, ``test_sn_cylindrical.py``, and
``test_sn_spherical.py``. Consolidating them into ``meshes.py`` is
deferred housekeeping (tracked in issue #77, "Reorganize tests/ by
model"); the module exists now so the eventual migration is a pure
search-and-replace against a stable import path.

Nexus integration
-----------------

sphinxcontrib-nexus ≥ 0.6.0 parses ``:math:`label``` roles in test
docstrings and writes a ``references`` edge from the containing
function/method node to the corresponding ``math:equation:*`` node. The
``verification_coverage`` and ``verification_audit`` MCP tools consume
these edges to build the test↔equation matrix.

**Requirements for the edge to form:**

1. The referenced label must exist as a Sphinx equation label (i.e.
   there is a ``.. math:: :label: collision-rate`` block in a
   ``docs/theory/**/*.rst`` page).
2. The test's containing file must be on Nexus's source path
   (``tests/`` is picked up automatically via
   ``nexus_test_patterns``).
3. The docstring must use the ``:math:\`label\``` form, *not* inline
   LaTeX source like ``:math:\`\Sigma_a\``` — the latter is correctly
   treated as inline math and produces no edge.

Rebuild Sphinx (``sphinx-build docs docs/_build/html``) to refresh the
graph. The MCP server reloads the database automatically on mtime
change, so a running agent session picks up the new edges without
restart.

Contributor checklist
---------------------

When adding a new test:

- [ ] Decide whether it is a **physics test** or a **foundation
  test**. Physics tests verify a ``:label:``\ -ed equation in
  ``docs/theory/**/*.rst`` and go on the L0..L3 ladder. Foundation
  tests verify a software invariant (data structure, numerical
  primitive, factory output) that has no theory label; they get
  ``@pytest.mark.foundation`` and **no** ``verifies(...)``. See
  :ref:`vv-foundation-tests` for the taxonomy and the anti-patterns.
- [ ] If it's a physics test, choose the right V&V rung. L0 is term
  verification against a hand calculation; L1 needs a *measured*
  convergence order; L2 is multi-group heterogeneous integration;
  L3 is experimental validation. 1-group tests are **degenerate**
  for transport — always demand ≥2G.
- [ ] Apply the level marker — ``@pytest.mark.l0`` / ... /
  ``@pytest.mark.foundation`` (or file-level ``pytestmark``).
  Don't rely on inheritance if the test isn't a thin wrapper around
  a single case.
- [ ] Physics tests: declare equation labels with
  ``@pytest.mark.verifies("label")`` and mirror them in the
  docstring as ``:math:`label``` so Nexus can link. If no theory
  label exists for what you're testing, the test is probably
  foundation — don't fabricate a label.
- [ ] If the test protects against a specific ERR-NNN or FM-NN, add
  ``@pytest.mark.catches("ERR-NNN", "FM-NN")`` and update
  ``.claude/skills/vv-principles/error_catalog.md`` to reference the new test by nodeid.
  The ``catches`` decorator is orthogonal to the level bucket — a
  foundation test can be the catcher for an ERR entry (ERR-020 is
  the canonical example).
- [ ] Run ``python -m tests._harness.audit`` and confirm your test
  appears in the expected level count. Run
  ``python -m tests._harness.audit --strict`` and confirm it still
  exits 0 (or the same exit code it had before your change, if the
  gate was already tripping on a pre-existing gap).
- [ ] If the test adds a new equation label to a theory page,
  rebuild Sphinx and confirm the Nexus graph has the edge via
  ``verification_coverage`` on that label. If the new label is
  definitional or points at code that doesn't exist yet, mark it
  ``:vv-status: documented`` per :ref:`vv-status-documented` —
  don't leave it as an orphan.
