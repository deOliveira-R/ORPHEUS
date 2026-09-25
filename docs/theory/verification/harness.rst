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
   :doc:`error_catalog`. No scattered assertions.

6. **Enforcement mode.** Every gate in ``tests/gates/`` carries a level
   tag — physics (``l0``..``l3``) or ``foundation``. The audit tool
   surfaces every untagged test and ``--strict`` exits non-zero on
   any gap, so new tests cannot slip in untagged. The "unmarked
   accumulates in its own bucket" stance from the initial migration
   is behind us; foundation finally gave the non-physics tests a
   home (see :ref:`vv-foundation-tests`).

7. **Type-error ratchet** (issue #226). The package carries a large
   pre-existing pyright error surface; until the per-module burn-down
   reaches zero, the enforceable invariant is monotonicity.
   ``tests/gates/test_pyright_ratchet.py`` (``foundation`` + ``slow``,
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

   Held by a withdrawal (N label(s) whose every carrier is a withdrawn
   test — neither covered nor orphan):
     <label>                                  N withdrawn (#NNN)

   Orphan equations (N of M testable theory labels have zero test
   coverage; K labels are .. vv-status: documented and excluded from
   the orphan gate):
     ...

   error_catalog.rst ERR coverage (N/M entries have a running catching test):
     DORMANT ERR-NNN (every catcher withdrawn: #NNN)

The output shape above is illustrative — every count drifts with each
commit, so run the CLI for current numbers (the auto-generated
verification matrix page carries the same data, refreshed on every
Sphinx build). The ``FD`` column counts
:ref:`foundation tests <vv-foundation-tests>`. Documented-only labels
are marked via the ``.. vv-status: <label> documented`` sentinel
described in :ref:`vv-status-documented`. The "Held by a
withdrawal" block and the ``DORMANT`` lines count tests marked
``@pytest.mark.withdrawn``, which are neither verifying nor catching
(:ref:`vv-withdrawn-generators`).

The tool runs ``pytest --collect-only`` under the hood so the
:data:`tests._harness.registry.TEST_REGISTRY` is populated, then
queries it. No test code is executed.

Flags:

``--json``
    Machine-readable output (full registry dump plus orphan /
    documented / phantom / ERR-coverage sets). Every relation is split
    into running and withdrawn carriers (``equation_coverage`` beside
    ``withdrawn_equation_coverage``, ``held_by_withdrawal``,
    ``err_coverage`` with ``running`` / ``withdrawn``,
    ``dormant_errors``, ``withdrawn_tests``), per
    :ref:`vv-withdrawn-generators`.
``--untagged``
    List only tests with ``level=None``. Should return an empty list
    under normal operation; non-empty output means new tests were
    added without a V&V tag (``l0``..``l3`` or ``foundation``).
``--gaps``
    List orphan equations (labels under ``docs/theory/**/*.rst`` with
    zero verifying tests, excluding documented-sentinel labels),
    phantom verifies-targets (tests naming a ``:label:`` that exists
    nowhere under ``docs/``), and ``ERR-NNN`` entries in
    :doc:`error_catalog` with no
    catching test. A label held by a withdrawal is not an orphan, and a
    dormant ``ERR-NNN`` entry (every catcher withdrawn) is not listed
    as missing a catcher; each is listed instead in its own section,
    "Held by a withdrawal" (label and issues) and "ERR entries held by a
    withdrawal" (entry and issues), so ``--gaps`` shows every gap the
    withdrawal is holding (:ref:`vv-withdrawn-generators`).
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
- There is usually **no** ``:label:`` in any ``docs/theory/**/*.rst``
  page that corresponds to what the test is checking, so most
  foundation tests carry no ``@pytest.mark.verifies(...)``. The
  level bucket is orthogonal to the equation link, though: a
  foundation gate that genuinely pins a labelled equation's content
  carries ``verifies`` **and stays foundation** — it is not promoted
  to L0..L3 (see the audit-reporting rules below).
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
from L0..L3. A ``verifies(...)`` mark is **orthogonal to the
level bucket**: the orphan gate counts marks from any collected
test, so a foundation gate that genuinely pins an equation's
content closes its orphan exactly as an L0–L3 test does (live
instances: ``tests/gates/sn/eigenvalue/test_keff_estimator_gate.py``,
``tests/gates/sn/operators/test_inverse_operator_equivalence.py``, the
reference-pillar ``V_αN`` gates). Most foundation tests carry no
``verifies(...)`` simply because they pin software invariants
that have no equation label; some algebra-of-record suites
additionally **opt out by declared design** (their docstring says
so) — those suites' labels carry a ``documented`` sentinel whose
rationale names the pinning gate. A label never carries both: a
live ``verifies`` target *is* tested, so a ``documented``
sentinel on it is a stale claim — drop the sentinel. Foundation
tests satisfy the ``--strict`` mode's untagged-tests gate — a
foundation marker is a valid tag.

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


.. _vv-withdrawn-generators:

Withdrawn reference generators — ``@pytest.mark.withdrawn``
-----------------------------------------------------------

A **reference generator** is a function under ``orpheus/derivations/``
that produces the value a test compares production against: an
eigenvalue, a flux shape, a matrix. A **withdrawal** is the
maintainer's ruling that a generator is not fit to be believed, so
that until the GitHub issue recording the ruling closes:

1. the generator is **not evidence**: a test that consumes it neither
   verifies an equation label nor catches a catalogued error, and the
   verification matrix and the error catalogue say so;
2. the generator **does not run**: a test that consumes it is skipped
   with the reason and the issue, and the generator itself refuses a
   call unless the run explicitly lifts the withdrawal.

Both consequences are needed. A skipped test whose ``verifies`` or
``catches`` marker still counted would report coverage that nothing
exercises; a counted-out test that still ran would spend its time on
a value nobody may cite. The generator and its tests are kept rather
than deleted, because the work that returns the generator to service
runs them (under the opt-in below) as its ladder.

The instance in the tree is `#506
<https://github.com/deOliveira-R/ORPHEUS/issues/506>`__: the Peierls
Nyström solver half of
``orpheus.derivations.continuous.peierls_nystrom`` (the volume kernel,
the boundary-closure operators, ``solve_peierls_1g`` and
``solve_peierls_mg``, the slab eigen-solve and the case builders the
lazy registry calls), withdrawn by the maintainer's ruling of
2026-09-24 because it is not research grade. The escape, transmission
and boundary primitives of the same package (``compute_P_*``,
``compute_G_*``, ``compute_T_*``), their two angular-assembly drivers,
the analytical identities of its ``reference`` module and the PS-1982
reference stay in service; :ref:`theory-peierls-nystrom` states which
half a symbol belongs to.

One value, read in two places
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A withdrawal is one value,
:class:`~orpheus.derivations.common.withdrawal.Withdrawal`, a frozen
pair ``(reason, issue)``. Its constructor refuses an empty reason and
an issue that is not a positive ``int`` (a ``bool`` or a numeric
string is refused too), and each refusal names the field. The value
is read in two places that must agree:

- **statically, by the test harness**, from the marker
  ``@pytest.mark.withdrawn(reason, issue=N)`` on every test that
  consumes the generator. ``tests/conftest.py`` parses the marker
  with :meth:`~orpheus.derivations.common.withdrawal.Withdrawal.from_mark`
  at collection time, skips the test, and records the
  ``Withdrawal`` on the test's registry entry
  (``TestMetadata.withdrawn``), where the audit reads it;
- **at run time, by the generator**, which is decorated with
  :func:`~orpheus.derivations.common.withdrawal.withdrawn_generator`
  and refuses the call (the lock, below).

The type lives in :mod:`orpheus.derivations.common.withdrawal`, in
production, because the lock is production code, and the conftest
imports the same type, so the marker's parser and the lock share one
definition. The module imports nothing outside the standard library
and nothing from ``pytest``: ``from_mark`` reads a marker only through
its ``args`` and ``kwargs`` fields.

For #506 the value is the package constant
``PEIERLS_NYSTROM_WITHDRAWAL`` (in
``orpheus/derivations/continuous/peierls_nystrom/__init__.py``): every
lock in the package is keyed on it. Its test-side marker is minted once
from the same value, as ``PEIERLS_NYSTROM_WITHDRAWN`` in
``tests/_harness/withdrawals.py``, and every #506 marker site uses that
one mark, so the reason sentence has one spelling. Gate M4 of
``tests/gates/test_withdrawal.py`` asserts that the mark parses to the
constant, that every #506 withdrawal in the collected tree carries the
constant's reason and issue, and (by an AST census of ``tests/``) that
no test spells ``pytest.mark.withdrawn(...)`` by hand.

The marker and where it may be placed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The spelling is one positional reason and one ``issue=`` keyword,
nothing else:

.. code-block:: python

   # tests/_harness/withdrawals.py mints the mark once:
   PEIERLS_NYSTROM_WITHDRAWN = pytest.mark.withdrawn(
       PEIERLS_NYSTROM_WITHDRAWAL.reason, issue=PEIERLS_NYSTROM_WITHDRAWAL.issue
   )

   # a test module uses it:
   from tests._harness.withdrawals import PEIERLS_NYSTROM_WITHDRAWN

   # on a module: every test in the file is withdrawn
   pytestmark = PEIERLS_NYSTROM_WITHDRAWN


   # on a class: every test method is withdrawn
   @PEIERLS_NYSTROM_WITHDRAWN
   class TestMGInputValidation:
       def test_rejects_mismatched_group_count(self): ...


   # on one test function
   @PEIERLS_NYSTROM_WITHDRAWN
   def test_slab_eigenvalue_converges(): ...

(The class and function names are illustrative.) pytest resolves
the marker as it resolves every marker, function first, then class,
then module ``pytestmark``, and the hook reads the closest one.

A withdrawal is a property of a test **function or class**, never of
one parametrised case. A ``withdrawn`` marker placed on a single
``pytest.param(..., marks=...)`` is refused at collection with a
``pytest.UsageError`` that names the test and says to withdraw the
function or class instead. So is a malformed marker: a missing
``issue=``, a second positional argument, any keyword other than
``issue``, or a value the constructor refuses. The refusal is a
collection error rather than a warning because a marker the hook
could not parse would leave the test running and counted. No
parametrised function in the #506 set mixes withdrawn and running
cases (``[M]`` 2026-09-25, the P0 specification's per-case
attribution, ``.claude/plans/reference_p0_spec.md`` §1.2), so a
per-case withdrawal is always a misplacement, never a need.

A withdrawn test is **skipped, never deselected**: the skip reason
reads ``withdrawn (#506): <reason>``, ``pytest -rs`` prints it, and the
default summary counts it under ``skipped``, so a run always shows how
many tests a withdrawal is holding. The complete #506 placement, one
node id per line, is ``tests/gates/withdrawal_506_placement.txt``;
gate M4 asserts that the tests withdrawn under #506 anywhere in the
collected tree are exactly that list, so adding or removing a marker
without updating the list is red, and so is the reverse.

The lock on the generator
~~~~~~~~~~~~~~~~~~~~~~~~~

The marker is a claim written on each consumer: "this test reaches a
withdrawn generator". Nothing in the marker checks that claim, and a
new test that calls the generator without the marker would run it and
count it. The check is placed on the generator itself, as a lockout
tag is placed on the breaker rather than on the equipment it feeds:
each withdrawn symbol is decorated with
:func:`~orpheus.derivations.common.withdrawal.withdrawn_generator`, and
a call raises
:class:`~orpheus.derivations.common.withdrawal.GeneratorWithdrawn`
unless the withdrawal's issue is lifted. The message names the
generator, the issue, the marker to add and the variable that lifts
it, for example::

   orpheus.derivations.continuous.peierls_nystrom.geometry.solve_peierls_1g
   is withdrawn (#506): Peierls Nyström reference is not research grade
   (maintainer ruling 2026-09-24). A test that reaches it carries
   @pytest.mark.withdrawn('Peierls Nyström reference is not research
   grade (maintainer ruling 2026-09-24)', issue=506); to run it anyway
   set ORPHEUS_RUN_WITHDRAWN=506.

The lock is what turns the placement into a checked property:

- **An unmarked test that reaches a withdrawn generator is red.** A
  marked test is skipped before it runs, so in a default run every
  test that does run is, by the lock, a test that reaches no
  withdrawn generator. Gate M5 runs every file that holds a #506
  withdrawal (and the one kept file beside them) under the canonical
  ``python -O -m pytest`` and asserts zero failures and exactly the
  listed skips. M5 is marked ``slow``, so ``-m "not slow"`` deselects
  it; every kept test in those files still enforces the property in
  its own run, since the lock fires wherever it is reached.
- **A reload cannot undo it.** The lock is part of the function
  object the module defines, so ``importlib.reload`` of the module
  re-executes the decorator and the reloaded symbol is locked again.
  A test-side tripwire that rebinds the symbols from a plugin is
  undone by a reload, and ``tests/gates/derivations/test_peierls_multigroup.py``
  reloads the Nyström ``cases`` module to re-read a routing flag
  (one of the two reasons the lock is in production, P0 specification
  §2.4).
- **A subprocess cannot escape it.** A worker process that imports the
  generator imports the lock with it, and the lift reaches the worker
  only through the environment. A pytest plugin never sees a child
  interpreter.

The decorator uses ``functools.wraps``, so the locked symbol keeps its
name, docstring, signature and ``__wrapped__``: ``inspect.signature``,
a registry keyed on names and a ``monkeypatch`` by attribute all see
the generator as before, and
:func:`~orpheus.derivations.common.withdrawal.withdrawal_of` reads the
lock back, so the set of locked symbols is enumerable. Gate M6
enumerates it twice: by an AST census of every
``@withdrawn_generator`` site under ``orpheus/``, and at run time over
the eight modules of the Nyström package, where each locked member
must be keyed on ``PEIERLS_NYSTROM_WITHDRAWAL``. Both must equal the
34 symbols of #506, so nothing else is locked (the primitives above
stay unlocked).

Why ``GeneratorWithdrawn`` is a ``BaseException``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:class:`~orpheus.derivations.common.withdrawal.GeneratorWithdrawn`
derives from :class:`BaseException`, not from :class:`Exception`, on
purpose. The refusal is a policy, not a failure of the computation, and
no fallback may absorb it. An ``except Exception`` around a generator
call (a root-finder loop that skips a failed bracket, a registry walk
that skips a broken producer, a test's ``pytest.raises(Exception)``)
would otherwise turn the lock into a silent fallback or a green test.
``KeyboardInterrupt`` and ``SystemExit`` are ``BaseException`` for the
same reason. The cost is that an ``except Exception`` cleanup block
does not run when the lock fires; a ``finally`` block does. pytest
reports the refusal as an ordinary test failure: its call wrapper
catches every ``BaseException`` and re-raises only its own ``Exit``
and ``KeyboardInterrupt``. The refusal pickles (``__reduce__``), so one
raised in a worker process reaches the parent as itself.

Why no ``xfail`` absorbs the refusal
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An ``@pytest.mark.xfail`` test reads ``xfailed`` when it raises,
whatever it raises (unless ``raises=`` narrows it), and that holds in
the setup phase too: a fixture that raises marks the xfail test
``x``. An unmarked xfail test that reaches a locked generator, in its
body or through a fixture, would therefore read ``x``, the same as the
expected failure the xfail documents, and the lock would enforce
nothing there. That population is real: ``[M]`` 2026-09-25, 12 of the
294 cases withdrawn under #506 were such xfails (recorded in the
hook's docstring). The ``pytest_runtest_makereport`` wrapper in
``tests/conftest.py`` closes it: when any phase (setup, call or
teardown) raised ``GeneratorWithdrawn`` and the report would read
xfailed, it rewrites the report to a failure whose message is the
lock's own (it names the issue and the marker to add). The wrapper is
tagged ``ELEGANCE-DEBT[guard] #506`` and retires with the lock. Gate M7
runs, in a child pytest, an unmarked xfail test whose body reaches a
stub lock, two whose function-scoped and module-scoped fixtures reach
it, and an ordinary xfail, and requires ``1 failed, 1 xfailed, 2
errors`` (a setup failure is counted as an error): the ordinary xfail
is the control leg that shows the wrapper rewrites only the lock's
refusal.

The opt-in: ``ORPHEUS_RUN_WITHDRAWN``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A withdrawal is lifted for one run by the environment variable
``ORPHEUS_RUN_WITHDRAWN``
(:data:`~orpheus.derivations.common.withdrawal.RUN_WITHDRAWN_VARIABLE`),
whose value is a comma-separated list of issue numbers, or ``all``:

.. code-block:: console

   ORPHEUS_RUN_WITHDRAWN=506 python -O -m pytest tests/gates/derivations/test_peierls_multigroup.py
   ORPHEUS_RUN_WITHDRAWN=506,512 python path/to/probe.py
   ORPHEUS_RUN_WITHDRAWN=all python -O -m pytest tests/gates/

With the issue lifted, the conftest hook does not skip the test and
the lock lets the call through. The variable is parsed by
:func:`~orpheus.derivations.common.withdrawal.lifted_withdrawals`, one
definition read by both; unset or empty lifts nothing, and an entry
that is not a positive integer is refused rather than read as "lift
nothing", because that would leave the caller believing the generator
had run. The conftest parses the variable once, in
``pytest_configure``, and a mistyped value (``ORPHEUS_RUN_WITHDRAWN=5o6``)
stops the session there with a ``pytest.UsageError`` naming the value.
The lock reads the variable at every call, so a variable set after
import takes effect, and a value that does not parse refuses the call
with a ``GeneratorWithdrawn`` naming the bad value rather than running
it. The ``Withdrawal`` value itself never reads the environment: each
reader asks ``withdrawal.issue in lifted_withdrawals()``.

It is an environment variable, and not a pytest command-line option,
for three reasons:

1. **A subprocess inherits it.** The worker process in
   ``tests/gates/cp/test_peierls_rank_n_protocol.py`` runs the
   generator in a child interpreter; a pytest option would need a
   second channel to reach it, and the environment is the channel the
   lock already reads.
2. **A script outside pytest can use it.** The improvement work under
   #506 runs the withdrawn generators from probes and scripts, which
   reach the lock and never reach a pytest option.
3. **Naming the issue scopes the lift.** ``ORPHEUS_RUN_WITHDRAWN=506``
   lifts #506 and nothing withdrawn under any other issue, so running
   one withdrawn family for its improvement work never silently runs
   another. A boolean flag could not say which withdrawal it lifts.

The variable **permits** a run and never changes a value: the
generator computes exactly what it computed before the withdrawal. A
lifted run's numbers are therefore the generator's numbers, still
withdrawn from evidence; the lift changes whether the code runs, not
whether its answer may be cited.

How the accounting counts a withdrawn test
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A withdrawn test is **neither verifying nor catching**. Its
``verifies`` and ``catches`` markers stay in the source, because they
describe what the test will claim again once the withdrawal lifts, but
the audit (``python -m tests._harness.audit``) splits every relation it
reads from the registry into **running** carriers (tests with no
withdrawal) and **withdrawn** carriers, and only running carriers
count as coverage:

- **Equation coverage** counts running carriers only. A label with at
  least one running carrier is covered by those; its withdrawn
  carriers are listed apart (JSON key ``withdrawn_equation_coverage``).
- **A label whose every carrier is withdrawn is held by a
  withdrawal**: neither covered nor orphan. The text report prints it
  under "Held by a withdrawal" with its withdrawn-carrier count and
  issue; the JSON carries it as ``held_by_withdrawal``
  (``{"carriers": N, "issues": [...]}``); the
  :doc:`matrix` renders it in its section "Claims held by a
  withdrawal". It is excluded from the orphan set, so ``--strict`` and
  ``--gaps`` do not report it as a gap: the gap is recorded, with its
  return criterion, in the issue the withdrawal names.
- **The phantom gate** (a ``verifies`` label that exists nowhere
  under ``docs/``) counts every carrier, withdrawn ones included: a
  dangling label is dangling whether or not its test runs.
- **Error catalogue coverage** is split the same way
  (``err_coverage``: ``running`` and ``withdrawn`` per entry). An entry
  whose every catcher is withdrawn is **dormant**, derived once in the
  audit and published as ``dormant_errors`` (entry to withdrawal
  issues), which the error index and the reconciler read: the text
  report prints
  ``DORMANT ERR-NNN (every catcher withdrawn: #506)``, the headline
  count reads "entries have a running catching test" and excludes it,
  and ``--gaps`` does not list it as missing a catcher.
- **Every withdrawn test** is listed in the JSON under
  ``withdrawn_tests`` with its reason and issue.

A dormant entry says so in its own body in :doc:`error_catalog`, on
the line

.. code-block:: rst

   **Status:** dormant — every catcher is withdrawn under #506

(several issues are listed ``#506, #512``). The reconciler's arm 7
(``test_a_dormant_entry_says_so_and_only_a_dormant_one`` in
``tests/gates/test_error_catalogue_reconciles.py``) holds the
catalogue to the audit's ``dormant_errors`` in both directions: a
dormant entry must carry the line, naming exactly its withdrawal
issues, and an entry that carries the line while it is not dormant is
a stale dormancy and is red. The arm reads pytest's
own marker resolution through the audit payload rather than parsing
decorators, so module, class and function placements are resolved
once, by pytest. Arm 2 of the same file (every ``catches`` claim names
a catalogued entry) is unchanged: a dormant entry's catchers are still
claims, and still name it.

The graph does not read the marker yet
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Nexus graph records ``verifies`` and ``catches`` edges from the
decorators whatever the test's other markers are, and it has no
notion of ``withdrawn`` (`deOliveira-R/sphinxcontrib-nexus#96
<https://github.com/deOliveira-R/sphinxcontrib-nexus/issues/96>`__).
So every graph answer counts a withdrawn test as a live carrier:
``nexus errors`` counts a dormant entry's catchers as coverage, and
``verification_coverage`` and ``verification_audit`` count a held
label as tested. Read a graph answer about a #506 label or entry
with that in mind; the audit and the matrix page are the sources that
know the split.

The generated error index (``.claude/skills/vv-principles/error_index.md``,
written by ``tools/verification/generate_error_index.py`` at the end
of every Sphinx build) reads the graph, so it carries an interim
correction: a **dormant** column and a "Dormant — every catcher is
withdrawn" section, read from the audit payload that the matrix
generator persisted earlier in the same build
(``docs/_generated/vv_audit.json``,
``tests._harness.audit.AUDIT_SNAPSHOT``), so the build collects the
suite once, not twice. The matrix generator deletes the previous
snapshot before it runs the audit and stamps the new one with the
audit's inputs (``HEAD``'s commit, whether the inputs differ from it
and a digest of the difference, untracked files included, plus the
collected count). The inputs are ``tests/``, ``orpheus/``,
``pyproject.toml`` and ``docs/theory/`` less the build's own generated
pages (the matrix page and the capability-matrix includes), which the
same build rewrites between writing the snapshot and reading it and
which feed nothing the audit reads. The index generator refuses, with
a message, a snapshot whose stamp is not the current inputs', rather
than using a stale one. Its catcher counts are still the graph's, so an
entry can read "5 catchers" and "dormant" on one row; the dormant
column is the one that says whether anything runs. The column retires
when the graph splits running from withdrawn carriers itself.

A transitional mechanism
~~~~~~~~~~~~~~~~~~~~~~~~

The lock is tagged ``ELEGANCE-DEBT[guard] #506`` in its docstring. It
is a run-time refusal standing where the reference machinery cannot
yet state "this reference is withdrawn" as a value. The
reference-solution architecture (#405, plan
``.claude/plans/reference_cache.md``) gives every reference solution a
certificate whose state is ``Valid``, ``Invalid`` (tripped
automatically by a failed check) or ``Withdrawn`` (a committed
declaration with a reason and an open issue; the generator never
runs). At its phase P4, when each family's generator returns a
certified reference, the ``withdrawn`` markers become ``Withdrawn``
certificate states and the decorator, the conftest hook and the
audit's split retire with them. Until then this section describes the
mechanism of record.

Nothing checks today that the issue a marker names is still open: a
withdrawal whose issue has closed keeps skipping its tests until
someone removes the markers and the locks together (the lock without
the markers reds every consumer, and the markers without the lock are
unenforced). Removing a withdrawal is one change that deletes the
locks, the markers, the placement list and the dormancy lines
together; ``[R]`` M4, M5, M6 and arm 7 between them red every partial
removal (markers without the list reds M4, a lock left behind reds its
now-unmarked consumers under M5, a changed lock set reds M6, a
dormancy line left behind reds arm 7).


Selecting tests at runtime
--------------------------

The standard pytest marker expressions apply:

.. code-block:: console

   pytest -m l0                       # only L0 term verification
   pytest -m "l1 and not slow"        # fast L1 checks
   pytest -m l2                       # integration (no gate carries l3)
   pytest -m foundation               # only foundation tests (software invariants)
   pytest -m "not foundation"         # only physics V&V
   pytest -m "l0 or foundation"       # L0 + foundation (fast; excludes eigenvalue runs)
   pytest -m "verifies and not slow"  # any test with an equation label
   pytest -m withdrawn -rs            # the withdrawn tests, each skipped with its reason

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
       references.py          # independent geometric references (mirror partners)
       mutation_batteries/    # hand-written, hazard-named mutation batteries

``xs.py`` re-exports the canonical cross-section helpers from
``orpheus.derivations.common.xs_library`` (``make_mixture``, ``get_mixture``,
``get_xs``, ``get_materials``, ``validate_all``) so tests can import
them from a single stable path. The Wigner-Seitz mesh builders
(``_ws_mesh``, ``_homogeneous_ws_mesh``, ``_two_region_ws_mesh``) are
defined privately in ``tests/gates/moc/test_verification.py``, and
``_build_homogeneous_mesh`` in ``tests/gates/moc/test_moc.py``; an empty
``meshes.py`` placeholder meant to collect them had no consumer and
retired with the move to ``tests/gates/`` (2026-09-22).

``tests/_harness/`` sits at the root of ``tests/``, beside
``tests/conftest.py`` and outside ``tests/gates/``, because every
regimen may import it, not only the pass/fail gates; the layout and
its reasons are :ref:`vv-test-suite-layout`.

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
2. The test's containing file must be on Nexus's source path. Every
   file under ``tests/`` is: the project-root pass analyses the whole
   tree (``.nexus/config.toml``; the ``nexus_*`` options that
   ``docs/conf.py`` once carried are retired).
3. The docstring must use the ``:math:\`label\``` form, *not* inline
   LaTeX source like ``:math:\`\Sigma_a\``` — the latter is correctly
   treated as inline math and produces no edge.

The graph does not read ``@pytest.mark.withdrawn``, so it counts a
withdrawn test's ``verifies`` and ``catches`` edges as live coverage;
see :ref:`vv-withdrawn-generators` for what that changes in a graph
answer and for the error index's interim dormant column.

Rebuild Sphinx (``sphinx-build docs docs/_build/html``) to refresh the
graph. The MCP server reloads the database automatically on mtime
change, so a running agent session picks up the new edges without
restart.

Contributor checklist
---------------------

When adding a new test:

- [ ] Decide its regimen first: a pass/fail check that must hold on
  every commit is a gate and goes in ``tests/gates/``, in the folder of
  the package it exercises; a timing or an experimental comparison is
  not a gate (:ref:`vv-test-suite-layout`, "Where a new case goes").
  The rest of this checklist is for gates.
- [ ] Decide whether it is a **physics test** or a **foundation
  test**. Physics tests verify a ``:label:``\ -ed equation in
  ``docs/theory/**/*.rst`` and go on the L0..L3 ladder. Foundation
  tests verify a software invariant (data structure, numerical
  primitive, factory output) that has no theory label; they get
  ``@pytest.mark.foundation``, and most carry no ``verifies(...)``
  because there is nothing to link. A foundation gate that genuinely
  pins a labelled equation carries ``verifies`` and stays foundation
  (the audit-reporting rules above). See
  :ref:`vv-foundation-tests` for the taxonomy and the anti-patterns.
- [ ] If it's a physics test, choose the right V&V rung. L0 is term
  verification against a hand calculation; L1 needs a *measured*
  convergence order; L2 is multi-group heterogeneous integration;
  L3 is experimental validation, whose cases belong in
  ``tests/validation/`` rather than among the gates. 1-group tests are **degenerate**
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
  :doc:`error_catalog` — add or extend the ``.. error-entry::`` block so it
  references the new test by nodeid.
  The ``catches`` decorator is orthogonal to the level bucket — a
  foundation test can be the catcher for an ERR entry (ERR-020 is
  the canonical example).
- [ ] If the test calls a withdrawn reference generator (for #506,
  the Peierls Nyström solver half), mark it with the withdrawal's
  minted mark (``@PEIERLS_NYSTROM_WITHDRAWN`` from
  ``tests/_harness/withdrawals.py`` for #506), on the function or class
  (never on a ``pytest.param``), and add its node id to the placement list
  (``tests/gates/withdrawal_506_placement.txt`` for #506). An unmarked
  test that reaches the generator is red with a message that names the
  marker (:ref:`vv-withdrawn-generators`).
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
