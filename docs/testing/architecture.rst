Test-Harness Architecture
=========================

.. contents::
   :local:
   :depth: 2

Motivation
----------

ORPHEUS operates a four-level V&V ladder. Each rung requires the rungs
below it:

.. list-table::
   :header-rows: 1
   :widths: 8 22 70

   * - Rung
     - Name
     - What it proves
   * - **L0**
     - Term Verification
     - Every individual term of the governing equation matches a hand
       calculation. Done with synthetic cross sections so the answer is
       closed-form. Catches sign flips, missing factors, wrong indices,
       and convention drift per term.
   * - **L1**
     - Equation Verification
     - The full operator reproduces an analytical eigenvalue, or the
       Method of Manufactured Solutions (MMS) residual converges at the
       asserted order. **Measured** convergence order, not just
       "passes a tolerance."
   * - **L2**
     - Integration
     - Multi-group heterogeneous problems self-converge or match an
       independent reference (cross-method consistency, Richardson
       extrapolation).
   * - **L3**
     - Validation
     - Comparison against experiment (ICSBEP, IRPhE).

1-group tests are **degenerate** for transport verification —
:math:`k = \nu\Sigma_f / \Sigma_a` regardless of flux shape, so angular
errors, normalization bugs, and convergence failures are invisible. The
ORPHEUS cylindrical diamond-difference recurrence bug (ERR-006 in
``tests/l0_error_catalog.md``) survived 20 one-group tests before a
multi-group run caught it. **Always demand ≥2G for any test claiming to
verify transport.**

Design principles
-----------------

1. **One source of truth per fact.** The V&V level a test belongs to,
   the equation labels it verifies, and the failure modes it catches
   are declared *once*, on the test itself, as pytest markers. The
   audit tool, the Sphinx verification-matrix page, and the Nexus
   knowledge graph all consume the same declaration.

2. **No new DSL.** The convention is vanilla ``pytest.mark.*``.
   :func:`tests._harness.verify` provides an ergonomic shortcut but
   teams that prefer raw markers keep full parity.

3. **Inherit from the reference case when possible.** Tests that pull
   analytical values via :func:`ref` inherit the V&V level from the
   underlying :class:`~orpheus.derivations._types.VerificationCase`.
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
   ``tests/l0_error_catalog.md``. No scattered assertions.

6. **Incremental migration.** The architecture is compatible with
   every existing test in the repository. Unmarked tests accumulate in
   the "unmarked" bucket of the audit report rather than breaking the
   build. Per-module migration PRs tag tests at their own pace.

Authoring a test
----------------

Preferred: the ``verify`` sugar layer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from tests._harness import verify


   @verify.l0(
       equations=["transport-cartesian"],
       catches=["FM-07", "ERR-003"],
   )
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

Effect of the decorator:

- every ``test_*`` method inside the class gets
  ``@pytest.mark.l0``, ``@pytest.mark.verifies("transport-cartesian")``,
  and ``@pytest.mark.catches("FM-07", "ERR-003")`` applied
- the docstring ``:math:`transport-cartesian``` role is picked up by
  sphinxcontrib-nexus and written as a graph edge from the test node
  to ``math:equation:transport-cartesian`` on the next ``sphinx-build``

Equivalent with raw markers
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import pytest


   @pytest.mark.l0
   @pytest.mark.verifies("transport-cartesian")
   @pytest.mark.catches("FM-07", "ERR-003")
   class TestSingleTrackAttenuation:
       ...

Both forms are first-class. ``verify`` is 30 lines of sugar; nothing is
hidden.

Parametrize over matching cases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Tests that exercise every matching :class:`VerificationCase` use
:func:`~tests._harness.verify.vv_cases`:

.. code-block:: python

   from tests._harness import vv_cases


   @vv_cases(level="L1", method="cp", geometry="slab")
   def test_cp_slab_eigenvalue(case):
       result = solve_cp_slab(case.materials, **case.geom_params)
       assert abs(result.k_inf - case.k_inf) < 1e-10

- Cases are pulled from
  :mod:`orpheus.derivations.reference_values` at collection time.
- The ``case`` parameter receives the full
  :class:`VerificationCase` instance.
- If every matched case shares the same ``vv_level``, the level marker
  is applied to the test automatically. Otherwise the test stays
  unmarked at the file level and individual parametrize IDs inherit
  their level from ``VerificationCase.vv_level`` via the conftest hook.

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

1. Explicit marker on the test itself (``@pytest.mark.lN``,
   file-level ``pytestmark``, or ``@verify.lN(...)`` which stamps an
   explicit marker).
2. Class name matching ``TestL<N>Foo`` (legacy convention, still
   honored).
3. Function name matching ``test_l<N>_*``.
4. Inherited from :class:`VerificationCase` via a parametrized
   ``case`` argument.
5. **Unmarked** — recorded in the registry with ``level=None`` so the
   audit tool can surface it.

Conflicts between different ``lN`` markers on the same test are
resolved deterministically (highest level wins) and a warning is
emitted so duplicates surface immediately.

The audit CLI
-------------

.. code-block:: console

   $ python -m tests._harness.audit
   ========================================================================
   ORPHEUS V&V Test Audit
   ========================================================================
   Total tests collected: 497

   By V&V level:
     L0            34   ( 6.8%)
     L1            15   ( 3.0%)
     L2             4   ( 0.8%)
     L3             0   ( 0.0%)
     unmarked     444   (89.3%)

   ...

Flags:

``--json``
    Machine-readable output.
``--untagged``
    List only tests with ``level=None``. Combine with ``grep`` to scope
    migration work per module.
``--gaps``
    List orphan equations (labels in ``docs/theory/*.rst`` with zero
    verifying tests) and ``ERR-NNN`` entries in
    ``tests/l0_error_catalog.md`` with no catching test.
``--strict``
    Exit 1 if any untagged tests or orphan equations are present.
    Used by CI in PR-final once migration is complete; **not** a CI
    gate in PR-1 so migration doesn't block merges.

    The ``--strict`` gate ignores any theory label that is marked
    :ref:`vv-status-documented` — those are deliberately excluded
    from the orphan set because they cannot or should not be paired
    with a test. A real gap (implemented-but-untested equation)
    still fires the gate.

.. _vv-foundation-tests:

Foundation tests — software invariants outside the L0..L3 ladder
-----------------------------------------------------------------

The L0..L3 ladder is defined around Cardinal Rule 4 — "Are we solving
the equations right?" Each rung assumes there is a **physics
equation** in a Sphinx theory page being verified: L0 checks a
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
   * - ``test_geometry::TestZoneSubdivision::test_equal_volume_*``
     - Every cell in an equal-volume zone has bit-identical volume
       by construction (the algebraic invariant that caught
       ERR-020). Not a physics claim — a round-trip
       correctness property of ``mesh1d_from_zones``.
   * - ``test_geometry::TestPWRFactories``
     - ``pwr_pin_equivalent(pitch)`` returns the correct cell
       radius from the Wigner-Seitz identity. Geometric primitive,
       not transport physics.
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
- There is **no** ``:label:`` in any ``docs/theory/*.rst`` page
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

Not every equation in ``docs/theory/*.rst`` can or should carry a
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

Rules:

- The ``vv-status:`` comment must appear in the same RST file as
  the ``:label:`` it refers to. Cross-file sentinels are not
  supported.
- The recognised status is **documented** only. Other words
  (``verified``, ``pending``, ...) are reserved for future use and
  are silently ignored today — they do not exclude the label from
  the orphan gate.
- If the label named in the sentinel does not actually exist in
  the RST (a typo), the sentinel is silently dropped and the real
  orphan continues to fire. Failing closed keeps typos visible.
- Do not use ``:vv-status: documented`` to paper over a genuine
  gap. "The test is hard to write" is not a justification;
  "the code does not exist yet" or "this is a definitional label"
  are. If in doubt, open an issue referencing the label.

The tool runs ``pytest --collect-only`` under the hood so the
:data:`tests._harness.registry.TEST_REGISTRY` is populated, then
queries it. No test code is executed.

Selecting tests at runtime
--------------------------

The standard pytest marker expressions apply:

.. code-block:: console

   pytest -m l0                    # only L0 term verification
   pytest -m "l1 and not slow"     # fast L1 checks
   pytest -m "l2 or l3"            # integration + validation
   pytest -m "verifies and not slow"  # any test with an equation label

Since ``verifies`` and ``catches`` are pytest marks with arguments,
``pytest -m "verifies"`` selects every test carrying any such mark.
Filtering by a specific label requires the audit tool (pytest's mark
expression language doesn't parse marker arguments).

``tests/_harness`` package layout
---------------------------------

.. code-block:: text

   tests/_harness/
       __init__.py      # re-exports verify, vv_cases, TEST_REGISTRY
       verify.py        # @verify.lN(...) decorators, vv_cases helper
       registry.py      # TestMetadata dataclass + TEST_REGISTRY dict
       audit.py         # python -m tests._harness.audit
       xs.py            # (stub) shared cross-section builders
       meshes.py        # (stub) shared mesh/geometry builders

``xs.py`` and ``meshes.py`` are stubs in PR-1. They expose the
eventual import surface so per-module migration PRs are pure
search-and-replace. The 200+ LOC of duplicated ``_make_pure_absorber_1g``,
``_ws_mesh``, etc. currently in ``test_moc_verification.py`` /
``test_cp_verification.py`` / ``test_sn_*.py`` are consolidated here
incrementally.

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
   ``docs/theory/*.rst`` page).
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

- [ ] Choose the right V&V rung. L0 is term verification against a
  hand calculation; L1 needs a *measured* convergence order;
  1-group is not enough for anything transport-related.
- [ ] Apply the level marker — ``@verify.lN(...)`` or raw pytest
  markers. Don't rely on inheritance if the test isn't a thin
  wrapper around a single case.
- [ ] Declare equation labels with ``@pytest.mark.verifies("label")``
  and mirror them in the docstring as ``:math:`label``` so Nexus can
  link.
- [ ] If the test protects against a specific ERR-NNN or FM-NN, add
  ``@pytest.mark.catches("ERR-NNN", "FM-NN")`` and update
  ``tests/l0_error_catalog.md`` to reference the new test by nodeid.
- [ ] Run ``python -m tests._harness.audit`` and confirm your test
  appears in the expected level count.
- [ ] If the test adds a new equation label, rebuild Sphinx and
  confirm the Nexus graph has the edge via
  ``verification_coverage`` on that label.
