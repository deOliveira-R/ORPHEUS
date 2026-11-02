.. _verification-cone-store:

=============================================================================
The cone store: the tests a working-tree change requires
=============================================================================

.. Machine header. Owner of: the cone of a test (the pins its verdict
   depends on) and the store that keeps one per test; why a cone is
   captured by fork with the test alone (the in-process cache leak), why a
   grandparent and one collecting parent per test module (ruling 38's
   refutation), the reading of a recording as a cone (roots, kept and
   dropped pins, the shadow presence pins, the cone-unknown reasons), the
   SQLite layout, the two-stage freshness, select's rule, the environment
   key, every read the capture silences and why, the positive controls and
   the mutation battery, the scope boundaries. Not owned here: the pin
   kinds, the recorder, validation, the kept reads, the fork guard and the
   refusals (reference_cache.rst, whose sections this page cites); the
   dependence DAG of the suite and its rests_on edges
   (principles.rst, vv-test-suite-layout); the V&V markers (harness.rst).
   Plan of record: .claude/plans/verification_tree_refinement.md, from
   "The cone store -- design" to its end (rulings 32 to 41, ruling 38 as
   corrected). Specification: scratch/verification_tree/cone_store/spec_b1.md;
   reviews qa_b1.md, elegance_b1.md, battery_b1.md, m1_diagnosis.md in the
   same folder. Written for increment B1 on branch feature/cone-store,
   uncommitted, 2026-11-01, over 6bc9a872.

The suite takes about an hour serially, and most edits change the verdict
of few tests. The cone store records, for each test, every node its verdict
depends on, and from those records lists the tests a working-tree change
requires. It is tooling, in ``tools/cone_store/`` (the user's ruling C6:
tooling under ``tools/``, data in SQLite under ``.cache/``), and it is the
second reader of the traced memo's machinery
(:ref:`verification-reference-cache`): it records with the memo's
``Recording``, names nodes with the memo's pin types, digests with the
memo's functions and judges freshness with the memo's ``Pin.stale`` and
file identities. It owns none of them (ruling 15 of the plan: one
change-detection layer).

.. code-block:: console

   python -O -m tools.cone_store capture [ids]   # record the cones of the named tests (default: all)
   python -O -m tools.cone_store select [ids]    # list the tests the working tree requires
   python -O -m tools.cone_store select --why    # ... with the reasons beside each id
   python -O -m tools.cone_store status          # what the store holds

Every command takes ``--repo`` (default the working directory), ``--store``
(default ``<repo>/.cache/cone_store/cones.sqlite``), ``--package`` (the
first-party package, default ``orpheus``) and ``--memo-root`` (the traced
memo's store, default the memo's own); ``capture`` takes ``--route
{auto,fresh}``. The ids are pytest node ids or paths, as pytest takes them.
The same three verbs are functions of the package (``tools.cone_store.capture``,
``select``, ``status``, and ``open_store`` for reading a store). The
command must run under the interpreter and ``-O`` level the suite runs
under, since both are in the environment key.

The full suite stays the merge gate. ``select`` is a local check and the
input of the census that places tests in the verification tree (ruling 22
of the plan); it replaces the full suite only after the validation run of
ruling 41 (sampled mutants, each run both in full and through the store,
every selection agreeing), which has not run.


Key facts
=========

* **A test's cone is the set of nodes whose change can change its
  verdict.** A node is one of the traced memo's pins: a def
  (``DefPin``), a module's import-time code (``ModulePin``, its skeleton),
  a data file's bytes (``DataPin``), what a path is (``PresencePin``), a
  directory's names (``ListingPin``), a memo entry read (``ChildPin``) or
  an installed distribution (``DistributionPin``). The interpreter, the
  platform and the rest of the environment are one key per cone, not pins.
* **A cone is captured with the test alone, by fork** (rulings 32 and 38).
  In one pytest process caches hide most of a cone: 127 of 177 in-suite
  cones of a 200-test sample missed defs the same test entered when run
  alone, a median of 35 defs each (`[M]` 2026-11-01). A child forked from a
  parent that has only collected gives the alone cone: 200 of 200 identical
  to a fresh interpreter's, at 0.29 s per test against 1.07 s.
* **Three levels of process.** A grandparent records from before the root
  conftest is imported; it forks one collecting parent per test module; each
  parent forks one child per test. A cache or ``sys.modules`` entry filled
  by one module's collection cannot hide code from another module's test,
  because the second module's parent never saw the first's collection.
  The grandparent loads every conftest a full run loads, so what they run
  (91 ``orpheus`` modules imported by the root conftest alone on the real
  tree, `[M]` 2026-11-01) is attributed to every test, an over-selection.
* **The store predicts full-run verdicts** (ruling 42). The full suite is
  the merge gate, and in a full run every conftest's module code and session
  hooks act on every test, so the capture reproduces that context: every
  conftest under the test paths is loaded before any test, and their set is
  in the environment key.
* **A cone is never smaller than the truth, by design; it may be larger.**
  Every rule that cannot decide drops toward selecting more: a dependency no
  pin can name makes the cone UNKNOWN and the test always selected; a test
  that did not pass at capture is always selected; a cone captured under
  another environment key is not trusted.
* **Freshness is the memo's two stages.** A pinned file whose identity on
  disk (device, inode, size, modification and change times) equals the
  identity stored at capture, settled then, is fresh unread; otherwise every
  pin on it is judged by its normalised-AST or byte digest, and a file found
  clean has its identity refreshed so the next ``select`` stops at stage one.
* **Over-selection is not bounded by any gate.** The ground-truth gate is
  containment (the store's cone holds the alone cone), which every superset
  passes. B2 (per-statement import-time pins, ruling 33) is what shrinks the
  grandparent's share.
* **Evidence.** 80 gate ids in 7 files (78 not ``slow``); a 35-arm mutation
  battery with a positive control, every arm reddening its target, and a
  14-arm round after review, every arm reddening its target and no
  declared-green row; the memo's own suites give the same outcomes captured
  and run plainly (M1, `[M]` 1342 s).


.. _verification-cone-store-cone:

What a cone is
==============

The cone of a test is the set of nodes whose change can change the test's
verdict. Its nodes are the traced memo's pins, so the store mints no node
type of its own; each is stale exactly when the memo's validation would call
it stale (:ref:`verification-reference-cache-manifest`):

.. list-table:: The pins a cone holds
   :header-rows: 1
   :widths: 20 42 38

   * - Pin
     - What it names
     - Stale when
   * - ``DefPin``
     - a first-party def, by its module path, qualname and ordinal (the
       index among same-named defs of the file), with its
       ``function_digest`` (normalised AST: docstrings and layout removed)
     - the def's normalised body changes, or the def or its file is gone
   * - ``ModulePin``
     - a first-party module's import-time code, by its
       ``skeleton_digest`` (constants, decorators, class attributes,
       imports: everything outside def bodies)
     - any edit to the module outside its def bodies
   * - ``DataPin``
     - a file the test read, by the digest of its bytes
     - its bytes change, or it appears or disappears
   * - ``PresencePin``
     - what a path the test probed is: a file, a directory, something else
       or nothing
     - that answer changes (a file created or deleted)
   * - ``ListingPin``
     - a directory the test listed, by the digest of its sorted names
       (``__pycache__`` left out)
     - a name is added, removed or renamed
   * - ``ChildPin``
     - a traced-memo entry the test read, by its key and payload digest
     - the entry no longer validates (its own manifest, checked recursively)
   * - ``DistributionPin``
     - an installed distribution whose code ran or whose file was read, by
       version and the digest of its ``RECORD``
     - another install of it is there

A row of the store is (test, environment key, pin), with the identity the
pinned file had at capture (ruling 39). The grain is the memo's at increment
B1: def bodies and module skeletons. An edit to a module's skeleton (a new
import, a constant) stales every test whose cone holds that module, whether
or not the test read the name the edit bound; per-statement pins (B2,
ruling 33) are the finer grain. The two directions are not symmetric. A
cone larger than the truth costs a test run; a cone smaller than the truth
lets a broken test go unselected. Every rule below is chosen so that a
case it cannot decide makes the cone larger.

A memo entry the test read is a ``ChildPin``, not the defs that generated
it: a warm hit runs no oracle code in the test's process. An edit to the
oracle stales the entry (its own manifest pins the oracle's defs), and the
stale entry stales its reader through the ``ChildPin``. ``select`` judges a
``ChildPin`` by the memo's validation, which reads and runs nothing and
never regenerates the entry (gate P2 asserts no process is started).


.. _verification-cone-store-capture:

How a cone is captured
======================

Why not in the suite's own process
----------------------------------

A capture inside one pytest run (a plugin recording which first-party code
each test enters) reads the cone of a test as the code the test RAN in that
process, and in a process that has run other tests much of that code has
already run and been cached. Measured `[M]` 2026-11-01 on sample N, the
200 tests of ``tests/gates/numerics/test_invariance.py``,
``test_measure.py`` and ``test_space_of_axes.py``, serial, ``-O``,
``OMP_NUM_THREADS=1`` (probes and logs in
``scratch/verification_tree/cone_store/``; the comparison covers 177 rows):

* 127 of 177 in-suite cones miss defs that the same test entered when run
  alone in a fresh process, a median of 35 defs per affected test and up to
  79; 0 of 177 in-suite cones hold a def the alone run did not enter;
* the missed defs are cache fills (``GeneratingMeasure.recurrence``,
  ``is_symmetric``, the ``Manifold`` methods);
* on the 87 non-slow tests of ``tests/gates/sn/eigenvalue``, two identical
  runs gave two cones for ``test_keff_slab.py::test_heterogeneous_absolute_keff``,
  according to whether ``StreamingCoefficientCache`` was built inside the
  test or found by ``__eq__``: an in-suite cone depends on run state, not
  only on test order;
* code a ``parametrize`` list runs (``gauss_legendre(8)``) runs at
  collection, outside every test's phases, and is attributed to no test.

The tree holds 189 decorated cache functions under ``orpheus/`` and 72
under ``tests/`` (`[M]` 2026-11-01, an AST census of decorators named
``cache``, ``lru_cache`` or ``cached_property``; the plan's 240 is a grep
for the same three names, a different predicate), beside hand-written
dictionary caches no census of decorators sees. Capturing with caches
switched off is refuted for the question "how is a cone captured": the
hand-written caches stay on, so its miss rate is unknown (ruling 32).

Why a fork, and why from a parent that has only collected
---------------------------------------------------------

A process that has done nothing but collect holds no cache a test filled.
A child forked from it runs one test and reports the recording: `[M]`
2026-11-01, sample N, 200 of 200 forked children passed, and 200 of 200
forked cones were identical to the cones of the same tests each run alone
in a fresh interpreter (setup plus call, first-party defs under
``orpheus/``). The fork costs 0.29 s per test (57.4 s for 200) where a
fresh interpreter costs 1.07 s (21.5 s for 20 tests of
``test_measure.py``) and the plain suite 0.03 s. Every level of the
capture is guarded by the memo's fork guard (below), since a fork after a
multi-threaded BLAS call crashed its child 20 of 20 times
(`[M]` 2026-10-16, :ref:`verification-reference-cache-fork-template`).

Why a grandparent and one collecting parent per test module
-----------------------------------------------------------

One collecting parent for the whole suite still leaks between test
modules, by two mechanisms, both in the parent:

* a cache filled while an EARLIER module was collected (a
  ``functools.cache`` filled by a ``parametrize`` list) is hit by a later
  module's test, which never enters the fill's defs;
* ``sys.modules``: a module's import-time code runs once, for its first
  importer, so its skeleton is attributed to whichever test module's
  collection imported it first.

`[M]` the test-architect, 2026-11-01, on a fixture tree: run in the suite
under one parent, the later test's cone lacks five nodes (the defs
``_build`` and ``table`` of the fill, and the module code of the three
modules its first importer imported); run alone it holds all five. On real
files (``test_measure.py``, ``test_invariance.py``,
``test_generating_measure.py``, 439 tests, 31 sampled and run alone), 31 of
31 sampled tests miss the module code of ``tests/gates/__init__.py`` and
``tests/gates/numerics/__init__.py``, run once at the first module's
collection. The 200 of 200 measurement above compared setup and call defs
under ``orpheus/`` only, so it could not see either mechanism.

So the capture forks one collecting parent PER TEST MODULE from a
grandparent that has imported only pytest and the conftests. Each module's
parent starts from the grandparent's state, never from another module's,
so neither mechanism can reach across modules. What the grandparent itself
imports or fills is attributed to every test, which over-selects and never
under-selects. Attributing the union of every module's collection cone to
every test is refuted for the same question: any skeleton edit anywhere
would select the whole suite.

.. dropdown:: First got wrong: "collection, recorded once per module in one parent, can only over-select"
   :color: muted

   Ruling 38 as first written read: code that runs at collection "is
   recorded in the parent once per test module and attributed to every
   test of that module, which can only over-select". The argument held for
   one module and failed across modules: the parent that collects module B
   after module A carries A's caches and A's ``sys.modules`` entries, so
   B's collection cone is missing whatever A's collection had already run.
   The fixture above is the positive control that refuted it, and gates C2
   and C5 (below) are its permanent form. The repair (grandparent, one
   collecting parent per module, one child per test) was ruled by the
   orchestrator under ruling 38's intent, the alone cone.

The processes as built
----------------------

``tools/cone_store/worker.py`` runs one plugin class per role, chosen once
from the role named on its command line:

1. **The grandparent** (``Grandparent``) starts the memo's ``Recording`` in
   a ``pytest_load_initial_conftests`` hook wrapper (``tryfirst``), before
   pytest imports the root conftest; ``pytest_configure`` would be too late,
   the conftest being imported by then. The recorder module is loaded by
   path (``tools/cone_store/_recorder.py``), registered under its dotted
   name ``orpheus.numerics._traced_memo_boot``: importing it by name would
   run ``orpheus/__init__`` and ``orpheus/numerics/__init__`` before the
   recorder started, and code that ran before the recorder is in nobody's
   cone. Pytest then imports its initial conftests, and the grandparent
   imports, recorded, every other conftest a full run would load
   (``tools/cone_store/conftests.py``: every ``conftest.py`` under the
   configured ``testpaths``, the rootdir when it names none, skipping hidden
   directories, ``__pycache__``, ``node_modules``, ``build``, ``dist`` and
   ``*.egg-info``), through pytest's own conftest loading. That walk itself
   runs unrecorded: it is the store's, and the environment key carries the
   set it finds. A conftest's module code and its session hooks
   (``pytest_configure``, ``pytest_collection_modifyitems``) act on every test
   of a full run whatever the test's directory: `[M]` qa 2026-11-01,
   ``tests/gates/sn/regression/conftest.py::pytest_configure`` sets the
   session's ``DriftWarning`` filter, and 5 test files outside that directory
   raise ``DriftWarning``. Once the conftests are imported, the grandparent forks one short
   child that writes how many pins the grandparent attributes to every test
   (``status`` prints it), then forks one collecting parent per test module,
   waiting for each.
2. **A collecting parent** (``Collector``) is that fork continuing as an
   ordinary pytest session over its module's arguments. It collects the
   module and, for each test, forks a child.
3. **A child** runs pytest's own protocol for the one test (setup, call,
   teardown), stops the recording, reads it as a cone
   (:ref:`verification-cone-store-reading`) and writes it, with its route
   and outcome, to the capture's results directory.

The recording is inherited across both forks and only grows, so a child's
cone is the grandparent's (the conftests' imports) plus its module's
collection plus its own run. Nothing re-arms the monitor: a code object the
recorder disabled after its first start is already in the set every later
child inherits.

**The fork guard.** Before every fork the process asks the memo's
:func:`~orpheus.numerics._traced_memo_boot.fork_blockers` with the platform
as it is now (:ref:`verification-reference-cache-fork-template`): a
platform outside ``FORK_MEASURED`` or more than one OS thread is a reason
not to fork. A collecting parent that may not fork runs each test in a FRESH
interpreter (``FreshChild``, route ``fresh``), which starts as the
grandparent does (no path argument, and the same load of every conftest) and
takes its node id once they are imported, so the two routes
record the same thing. A grandparent that may not fork spawns a
``FreshCollector`` per module, which records nothing itself and runs every
test fresh; ``--route fresh`` does the same for every module. Every fork and
spawn of the capture is made inside the memo's ``spawning()``, so it is not
recorded as the test's own process start.

**What a child reports.** Its route (``fork`` or ``fresh``) and its outcome:
``passed``, ``failed`` or ``skipped`` from pytest's reports; ``crashed``
when a fork child or a fresh interpreter produced no result (the cone is
then unknown, with the wait status as its reason); and, written by the
driver, ``uncaptured`` for a collected test no worker reported and
``collection-error`` for a module that failed to collect.

**The bytecode cache.** Every worker sets ``sys.pycache_prefix`` to
``<repo>/.cache/cone_store/pycache`` inside its own process and never
exports ``PYTHONPYCACHEPREFIX``. The prefix keeps the capture from reusing
bytecode compiled in another checkout (a code object is attributed by its
file name, and pytest's assertion-rewriting cache does not repair the name
on load). The variable is not exported because every interpreter a test
starts would inherit it: `[M]` 2026-11-01 (the test-architect's M1
diagnosis), with it exported, ``test_traced_memo_source_index.py::test_b17``
flipped from pass to fail under capture, since it deletes
``orpheus/numerics/__pycache__`` to force a child to compile fresh and that
directory no longer existed. A capture must not change the system under
test.

**Collection.** ``select`` and ``capture`` both start with one
``--collect-only --continue-on-collection-errors`` subprocess (role
``list``), which writes the node ids pytest would run, in its order, and the
modules whose collection failed. That collection is ``select``'s cost floor.


.. _verification-cone-store-reading:

Reading a recording as a cone
=============================

``cone_of`` (``tools/cone_store/cone.py``) turns a stopped recording into a
cone in three steps.

**Identities first.** For every first-party code file and every
repository data file the recording names, the file's identity is read by
the memo's ``Witness``, with ``reference_ns`` set to the capture's start,
BEFORE any pin reads the file's bytes (the memo's stat-before-read order,
:ref:`verification-reference-cache-kept-reads`). A file last written less
than ``TRUSTED_AGE_NS`` (2 s) before the capture began, or during it, is
stored unsettled.

**The pins.** :func:`~orpheus.numerics.traced_memo.collect_pins` returns
the pins and every refusal as values
(:ref:`verification-reference-cache-refusals`). The cone keeps:

* ``DefPin`` and ``ModulePin`` of first-party code, a file under
  ``<repo>/<package>/`` or ``<repo>/tests/``; a test's own body and its
  helpers under ``tests/`` are first-party, so an edit to a test selects it
  at def grain;
* ``DataPin``, ``PresencePin`` and ``ListingPin`` of paths under the
  repository, except under the bytecode cache directory and the directories
  above it that the capture itself creates;
* ``ChildPin`` and ``DistributionPin`` as they come.

It drops, by decision, the interpreter, environment, working-directory and
memo-source pins: the environment key covers the first two, and the memo's
own sources are first-party files whose defs and skeletons are already
pinned. All paths are stored relative to the repository, so a store
survives a move of the checkout; ``tools/cone_store/pins.py`` is the one
function that changes a pin's path, and the stored row of a pin is the
memo's own (``Manifest.to_json``: ``[kind, *fields]``, decoded by
``Manifest.from_json``).

**Shadow presence pins.** The import system's own stats and listings are
exempt from the recording (they are how a module is found;
:ref:`verification-reference-cache-manifest`, "The import system is
exempt"). So a new module or package that would shadow an imported one is
seen by no recorded pin. For each first-party module that ran, ``cone_of``
adds a ``PresencePin`` for every name that would be found before it: in its
own directory, a package directory or an extension module of the same name;
for a top-level name, the same names in the other ``sys.path`` entries
under the repository. They are presence pins and not a listing of the
directory, because a listing would move whenever an unrelated file is
added there. Gate S12 (a new package ``pkg/consts/`` shadowing
``pkg/consts.py`` selects the importers) is their witness.

**Cone unknown.** Each refusal is read in one exhaustive ``match`` on its
kind; a kind it does not know raises ``TypeError``, so a new kind added to
the memo cannot be ignored silently:

.. list-table:: When a cone is unknown
   :header-rows: 1
   :widths: 30 70

   * - Cause
     - Rule
   * - ``UnadmittedProcess``
     - the test (or its module's collection, or the grandparent) started a
       process; that process's reads are unrecorded (ruling 40). Reason:
       ``process started: <program>``
   * - ``ChangedByRun``
     - a repository file the run read and then changed; a path outside the
       repository (a temporary file the test made and read back) is dropped
   * - ``NoRecordedBytes``, ``SourceGone``, ``NoDefAtLine``
     - when the file is first-party; dropped otherwise
   * - ``AmbiguousOwner``, ``UnownedSitePackage``, ``UninstalledDistribution``
     - always
   * - changed during the capture
     - a pinned file that is not settled against the capture's start and
       whose identity differs from the one it had when the capture began:
       the pins read the file's bytes after the test ran, so they may pin
       bytes that never ran
   * - a crash
     - a fork child or fresh interpreter that wrote no result

Every reason is stored as text beside the test, never as a bare flag, so
``status`` and ``select --why`` say what to look at.


.. _verification-cone-store-store:

The store
=========

``tools/cone_store/store.py``, SQLite at ``.cache/cone_store/cones.sqlite``
(``.cache/`` is ignored by git):

.. list-table::
   :header-rows: 1
   :widths: 14 86

   * - Table
     - Contents
   * - ``tests``
     - one row per test id: the environment key it was captured under, its
       route, its outcome and its cone-unknown reasons (a JSON list)
   * - ``pins``
     - each distinct pin once, as the memo's row, with the file it judges
       when it is a content pin (``DefPin``, ``ModulePin``, ``DataPin``)
   * - ``files``
     - one row per judged file: its identity at capture and whether that
       identity is trusted
   * - ``cones``
     - (test id, pin id): the cone of each test
   * - ``meta``
     - the grandparent's pin count and the number of tests the last capture
       pruned
   * - ``rows`` (a view)
     - the flat reading, one row per test and pin, with the file's identity

**A file's identity is a fact about the file**, stored once and shared by
every pin on it. A capture that reads a file whose stored identity differs,
or was not settled, stores it untrusted, so the next ``select`` judges every
pin on that file (stage two). **A capture replaces a test's rows whole.** It
forgets tests that no longer exist: those of a re-captured module that the
module no longer holds, and, after a capture of everything, those of a module
that is gone. A capture of named node ids prunes nothing else, and that is
safe: ``select`` iterates collected ids only, and a re-added test's old row
is judged by its pins' digests against the current tree.


.. _verification-cone-store-freshness:

Freshness in two stages
=======================

``select`` judges each pin of each candidate test as the memo judges an
entry (:ref:`verification-reference-cache-kept-reads`), file by file:

1. **Stage one, the file's identity.** The pinned file is stat'ed (before
   anything reads it). When its identity equals the stored one and the
   stored one is trusted, every pin on the file is fresh, and the file is
   not opened.
2. **Stage two, the pins.** Otherwise every pin on the file is judged by the
   memo's own ``Pin.stale`` against one ``_Checkout`` (``traced_memo``'s
   reader of first-party sources), so a def pin compares the normalised-AST
   digest and a docstring, comment or blank-line edit stays fresh. When
   every pin on the file is clean, the file's stored identity is refreshed
   with the identity just stat'ed and whether it is settled, so the next
   ``select`` stops at stage one; a file with a stale pin is judged again
   every time.

Pins that judge no file (``PresencePin``, ``ListingPin``, ``ChildPin``,
``DistributionPin``) are judged directly. The memo resolves a ``DefPin`` or
``ModulePin`` path through ``sys.path``, so ``select`` puts the repository
first on ``sys.path`` while it judges.

**Why the identity is trusted only when settled.** A second write within
one timestamp tick can leave a file's size and both times as the first
write left them (git's "racily clean" case). An identity is trusted only
when the file was last written more than 2 s before the moment it is
judged against: for a capture, the capture's START, not the moment of the
stat. A capture can run for an hour, and a file edited at its first minute
would be settled against a stat taken at its last.

Gates F1 to F6 pin this: F1, no file changed, ``select`` opens no source and
builds no source index, and after one file is touched exactly that file is
opened; F2, a file written just before the capture is read again at stage
two; F3, an operator swap with the size and modification time restored by
``os.utime`` is seen (the change time moves on every write and no call
restores it); F4, a comment edit passes stage one and stops at stage two; F5,
a second ``select`` in a new process after F4 builds no source index (the
refresh); F6, a def moved down the file stays fresh.


.. _verification-cone-store-select:

What ``select`` lists
=====================

``select`` collects the current tests, then lists, in pytest's order:

1. every module whose collection failed, by its path ("collection error
   in <module>");
2. every collected test the store does not hold ("new: not in the store");
3. every test captured under another environment key;
4. every test whose cone is unknown, with each reason;
5. every test that did not pass at capture: its cone stops at the failure,
   so it is partial ("did not pass at capture (failed)");
6. every other test with at least one stale pin, with the reason each pin
   gives (the def or module gone, the digest moved, the file's bytes
   changed, the child no longer valid, the path's presence changed).

A test that is no longer collected is not listed. A test parametrised anew
is new under its new ids. A decorator edit is a skeleton edit. With no edit
after a capture, ``select`` lists exactly the unknown tests and the tests
that did not pass (gate S1). ``status`` prints the store's counts: tests by
route, cone-unknown, collection errors, not passed, under another
environment key, distinct pins and rows, the grandparent's pin count, the
last capture's pruned count, and each unknown test with its reasons.


.. _verification-cone-store-environment:

The environment key
===================

A cone is valid only under the environment it was captured in (ruling 34),
so each test row carries the digest of:

* the interpreter, its version string and bytecode tag
  (``traced_memo._python_identity``);
* the platform: the operating system, the machine, the ``-O`` level and
  numpy's BLAS (:func:`~orpheus.numerics.traced_memo.platform_tag`);
* every installed distribution with its version;
* the environment variables a cone can follow: the memo's own list of the
  variables it fixes for a generation (``PYTHONHASHSEED`` and the five
  thread-count variables ``OMP_NUM_THREADS``, ``OPENBLAS_NUM_THREADS``,
  ``MKL_NUM_THREADS``, ``VECLIB_MAXIMUM_THREADS``,
  ``NUMEXPR_NUM_THREADS``; one list, the memo's ``_FIXED_ENVIRONMENT``,
  read and never re-spelled), plus ``PYTEST_ADDOPTS`` (it changes what
  pytest runs), ``PYTEST_PLUGINS`` and ``PYTEST_DISABLE_PLUGIN_AUTOLOAD``
  (they change which plugins load), ``PYTHONWARNINGS`` (it can turn a warning
  into a failure) and ``MKL_CBWR`` (it changes how MKL sums);
* every installed ``pytest11`` entry point, as (name, value, distribution):
  pytest loads these before any recorder starts, and a distribution can gain
  one under an unchanged version (gate S21);
* the set of conftests a full run loads (``tools/cone_store/conftests.py``,
  the same function at capture and at ``select``), so a conftest that appears
  or disappears anywhere selects every test (gates S13, S19);
* the bytes of the repository's pytest configuration files,
  ``pyproject.toml``, ``pytest.ini``, ``setup.cfg`` and ``tox.ini``: pytest
  reads them before any recorder starts, so no pin could see an edit to a
  ``filterwarnings`` or ``addopts`` setting (qa's review, U4; gate S14).

``select`` lists every test whose row carries another key. A wider key only
over-selects: any ``pip install`` of an unrelated package moves the
distribution table and re-selects the suite. Whether the per-test
``DistributionPin`` rows could replace the whole table is open (the
elegance review's question C2); ruling 34 names "the distributions" and the
key holds them all.


.. _verification-cone-store-silenced:

What the capture does not record, and why
=========================================

A recording offers a test everything its process does, and some of it is
pytest's work or the store's own, not the test's dependency. Recorded, each
of these made a whole-suite selection out of an edit that changes no
verdict. Each silenced read runs under the memo's ``unrecorded()``, the flag
the recorder's audit hook and its ``os.stat`` probe both honour, so the
same read made by a TEST (a census reading a source, a docstring gate) is
still recorded.

.. list-table:: Silenced reads
   :header-rows: 1
   :widths: 30 70

   * - Read
     - Why it is not the test's
   * - pytest's assertion rewriter reading each test module
       (``_pytest.assertion.rewrite._rewrite_test``)
     - it compiles the module the test is defined in; the module's code is
       pinned where it runs, by def and skeleton. Recorded, the whole-file
       read pinned the test module's bytes, so a comment edit selected it.
       A test that reads a source as text is still recorded: gate S11, a
       test that reads ``pkg/user.py`` and is selected by a docstring edit
       to it
   * - collection's directory scans (``scandir`` in ``_pytest.main`` and
       ``_pytest.python``) and its virtual-environment probe of each
       directory entry (``_pytest.main._in_venv``)
     - how pytest finds the test files; a new test file changes no existing
       test. `[M]` 2026-11-01 (battery, S8): recorded, ``ListingPin('tests')``
       and pytest's ``_in_venv`` probes of every entry
       (``tests/_helpers.py/pyvenv.cfg``) stood in every cone, so deleting
       any file under ``tests/`` selected every test
   * - pytest's own implementations of the collection walk's
       ``pytest_ignore_collect`` hook (an ``is_dir`` probe of every entry of
       every directory it descends, even in a session over one module)
     - the same: only the implementations whose plugin module is
       ``_pytest.*`` run unrecorded, so a sibling file's presence is no
       test's dependency, while a conftest's or a plugin's own
       ``pytest_ignore_collect`` stays recorded (gate S20: a user hook that
       fills a cache from a file is selected by an edit to that file)
   * - pytest's initial conftest discovery
       (``PytestPluginManager._try_load_conftest``, which globs
       ``<testpath>/test*`` and asks ``is_dir()`` of each entry)
     - in the grandparent, so it was attributed to every test: deleting any
       test file selected the suite. The loading of each directory's
       conftests inside it (``_loadconftestmodules``: the probe of
       ``conftest.py`` and the import) is re-enabled by a recorded wrapper;
       a NEW ``conftest.py`` anywhere is seen through the environment key's
       conftest set (gates S13, S19)
   * - the store's own walk for every conftest a full run loads
     - its result is the environment key's conftest set, not a test's
       dependency
   * - the memo's digest of its own source (``_memo_digest``, run when
       ``traced_memo`` is imported, inside the capture)
     - the memo's own read, under ``unrecorded()``; offered, the memo's whole
       bytes would be a data dependency of every test
   * - distribution discovery
       (``importlib.metadata`` and, when installed, the
       ``importlib_metadata`` backport: ``FastPath.children``)
     - it lists every ``sys.path`` entry for metadata directories; which
       distributions are installed is the environment key's fact.
       Recorded, the listing of the repository root and of ``tests/`` made
       every added or deleted file there select every test
   * - the capture's own forks and spawns
     - made inside the memo's ``spawning()``, so they are not the test's
       process starts (gate U7: no non-spawning test is unknown); likewise
       the memo's own generations (gate U8)
   * - ``trace_call``'s manifest re-read
     - the memo's own read of the traced sources, made under
       ``unrecorded()`` inside ``trace_call``
       (:ref:`verification-reference-cache-unseen`); gate M2 and the
       memo's nested-recording gate

What the cone drops after recording (the interpreter, environment,
working-directory and memo-source pins, the bytecode cache directory) is
listed in :ref:`verification-cone-store-reading`. The memo's own exemptions
(the import system's reads, by the immediate caller; the memo's store reads)
hold in the capture as in a generation.


.. _verification-cone-store-controls:

The positive controls and the mutation battery
==============================================

The plan wrote four positive controls before the build (X1 of the
instrument doctrine); each is a gate on a fixture tree that
``tests/gates/tools/_cone_store_fixture.py`` writes under ``tmp_path`` (a
package ``pkg/`` with a constant module, a cached fill, a module with twin
defs, a traced-memo oracle, and test modules each exercising one mechanism,
with variants built by parameters of the one builder; the package is not
named ``orpheus`` so the editable install cannot shadow it):

1. **A constant edit selects its readers** (P1, with C3): at B1's module
   grain, the edit to ``TOL`` in ``pkg/consts.py`` selects exactly the
   tests whose cones hold ``pkg/consts.py``'s skeleton, a superset of the
   readers of ``TOL`` that is exact on the fixture. It becomes exact by
   construction in B2.
2. **An oracle edit selects its validation test and its readers** (P2): the
   edit to the oracle's body selects the test that runs the oracle (by its
   ``DefPin``) and the test that reads the warm entry (by its ``ChildPin``),
   and ``select`` starts no process.
3. **A prose or layout edit selects nothing** (P3, four ids): a def
   docstring, a module docstring, a comment line, a blank line between two
   defs.
4. **Ground truth** (T1 to T4): every fixture test's stored cone contains
   its cone captured alone by ``--route fresh``, with 0 misses over every pin
   kind (T1); the same on a real-suite sample chosen by mechanism, 14 tests
   in five files plus one spawner asserted unknown (T2); the same on sample
   N's 200 tests (T3, ``slow``); the capture's and the alone run's outcomes
   agree (T4). The alone capture shares the recorder with the store, so the
   comparison is independent on the route axis only (X4); the recorder's own
   ground is R1 to R5, whose expected sets are written by hand from the
   fixture's source.

The gates, ``tests/gates/tools/test_cone_store_*.py``, every one
``foundation`` (`[M]` 2026-11-02, ``pytest --collect-only``: 80 ids, 78 not
``slow``):

.. list-table::
   :header-rows: 1
   :widths: 26 14 60

   * - File
     - Rows
     - Claims
   * - ``test_cone_store_recorder.py``
     - R1 to R9
     - exact pin sets against hand-written expectations; every digest is
       the memo's; the store defines no digest (R7), starts no second
       recorder (R8) and asks ``fork_blockers`` rather than re-spelling the
       guard (R9)
   * - ``test_cone_store_collection.py``
     - C1 to C6, H1
     - a collection fill is attributed to its module; another module's fill
       cannot hide it (C2); the grandparent's code reaches every test (C4);
       a package ``__init__`` reaches every test under it (C5); the rows do
       not depend on module order (C6); the in-suite monitor's blindness to
       a fill is asserted, not assumed (H1)
   * - ``test_cone_store_selection.py``
     - S1 to S22, P1 to P3
     - ``select``'s rule, the three controls, the environment key per
       component (S5, S6), the shadow pins (S12), a new conftest selects
       every test (S13, S19), the ini bytes (S14), a collection error (S15),
       a read-then-rewritten file (S16), an edit during the capture (S17),
       an off-path conftest's session hook (S18), a user
       ``pytest_ignore_collect`` stays recorded (S20), a ``pytest11`` entry
       point under an unchanged version (S21), ``PYTEST_PLUGINS`` (S22)
   * - ``test_cone_store_freshness.py``
     - F1 to F6
     - the two stages, by a spy on the files opened and the source indexes
       built (a claim about the path is gated on the path)
   * - ``test_cone_store_routes.py``
     - G1 to G4, U1 to U9
     - the fork guard at the grandparent and at a module, the unmeasured
       platform, fork equal to fresh (G4); each process API marks a test
       unknown (U1 to U5: ``subprocess.run``, ``os.system``, a spawn-context
       ``multiprocessing`` process, ``os.fork``, ``os.posix_spawn``), at
       collection too (U6), and the store's and memo's own processes do not
       (U7, U8)
   * - ``test_cone_store_ground_truth.py``
     - T1 to T4, CU-OVER
     - control 4; CU-OVER records the grandparent's pin count, a figure
       and not a bound
   * - ``test_cone_store_memo_coexistence.py``
     - M1 to M3
     - the memo's own suites give the same outcomes captured and plainly
       (M1, ``slow``); a ``trace_call`` nested in a capture is seen by it
       and not polluted by it (M2); the refusals are values (M3)

**The battery** (`[M]` the test-architect, 2026-11-01,
``scratch/verification_tree/cone_store/battery_b1.md``): one textual
mutation per arm, installed in a worktree copy with a check that its target
text occurs exactly once, the arm's rows run under ``python -O -m pytest -m
"not slow"`` with the worktree on ``PYTHONPATH`` and each file restored from
a pristine copy (``diff -rq`` empty after the last arm). The positive control
COLL-DROP (each child's recording cleared before its test) reddened six ids
over five rows, so the harness was live. 35 arms ran and 0 reddened no row.
Among them: COLL-SHARED (one parent for every module) reddens R1, C3 and G4;
GP-DROP (the grandparent's record not attributed) reddens C4, C5, R1 and T1;
REC-LATE (recording started at ``pytest_configure``) reddens C4 and R1;
ROOT-T (``tests/`` out of the roots) reddens R1, R2 and S8; each guard bypass
reddens its own G row; FRESH-IDENT (freshness by file identity alone)
reddens P3, F4 and F5 with over-selection messages only; IDENT-WEAK (identity
as size and modification time) reddens F3; ORD (the ordinal dropped from a
def's key) reddens R5. After review, 14 more arms (shadow pins dropped,
presence pins dropped, the ini bytes out of the key, collection errors
ignored or unrecorded, a repository read-then-write not unknown, the
changed-during-capture check skipped, the ``.py`` data-pin drop restored, an
own digest, an own guard, and four arms on the memo copy: ``origin()``
outside ``pinned``, one refusal kept, ``trace_call``'s manifest recorded,
``_LISTING_IGNORED`` emptied) each reddened its target row and no
declared-green row.

**Cells that stay green by design**, each with its structural reason:

* T1 under COLL-DROP and REC-LATE: the alone capture runs through the same
  recorder hook, so the mutation lands on both sides (the shared-recorder
  axis, declared); R1 and C1 to C5 are the catchers.
* COLL-SHARED leaves C2, C5, C6 and T1 green: one recording accumulating
  from the grandparent attributes the UNION of all collections to every
  test, which over-selects; the under-selection of ruling 38 as first
  written needed a reset per module, which is COLL-DROP plus sharing.
* ROOT-T leaves S2 and S10 green: with ``tests/`` outside the roots, a test
  file's whole-file data pin is kept and selects on any edit, so the store
  degrades to file grain for test files instead of failing.
* NOFORK (each test run in the collecting parent) leaves H1 green, because
  the recording is cumulative, and reddens G4, T1 and the fixture's T4.
* Not run, each for a structural reason: MON-NOREARM (the package holds 0
  ``restart_events`` calls: nothing to mutate, since the recording only
  accumulates); UNK-MEMO (the memo's ``spawning()`` flag is the same check
  that exempts the store's own forks, so the arm is UNK-SELF and UNK-MEMO at
  once); READS-REGEN (the store never calls the memo's ``__call__``; P2's
  route leg guards it); FRESH-FILE and KEY-LINE (not one textual edit; S7 and
  F6 are their catchers).

**Runs** (`[M]` 2026-11-01, before the full-run repairs of ruling 42):
``python -O -m pytest -m "not slow" tests/gates/tools/`` (the cone store's
gates and the other tools gates): 176 passed, 0 failed, 606 s; M1 alone
(``slow``): 1 passed, 1342 s. The directory now collects 181 non-slow ids
with S18 to S22; their run is the orchestrator's, after this page.
S18 to S22 and the re-posed S13 are written from qa's second pass
(``scratch/verification_tree/cone_store/qa_b1_delta.md``), each with the arm
named in its docstring as its first red.


.. _verification-cone-store-limits:

Scope boundaries
================

Each is tagged ``SCOPE-BOUNDARY[guard]`` in ``tools/cone_store/__init__.py``
with the machinery that would lift it:

* **Environment reads beyond the key are invisible.** No audit event fires
  on a read of ``os.environ``, so a test whose verdict follows a variable
  outside the key (``ORPHEUS_RUN_WITHDRAWN``, say) is not re-selected when
  it changes. The machinery is a recorder of environment reads (#687).
* **A pin reads its file's bytes after the test ran.** A capture child has
  no pinned-source loader, so ``collect_pins`` reads each first-party source
  after the run, as ``trace_call`` does. A file written within 2 s before
  the capture began, or during it, is stored unsettled, and a first-party
  file whose identity at the end of a test differs from its identity when
  the capture began makes that test cone-unknown (gate S17). A file edited
  and restored, or edited after its test and before the test is judged,
  within one timestamp tick of the stat, is not seen. The machinery is the
  memo's pinned-source loader extended to pytest's assertion-rewriting
  loader (#688): a B1 scope boundary and a B2 candidate, ruled by the
  orchestrator on 2026-11-01.
* **Over-selection is unbounded.** Whatever the grandparent imports is
  attributed to every test: the root conftest imports 91 of the 403 ``.py``
  files under ``orpheus/`` (`[M]` 2026-11-01, a fresh ``python -O`` importing
  ``tests.conftest``, counting ``orpheus`` modules in ``sys.modules``), so
  an edit to the skeleton of any of those modules selects the whole suite. A
  new file selects every test that probed or listed its directory. The
  machinery is B2's per-statement pins (#689).
* **A process a test starts is unrecorded.** Such a test is cone-unknown and
  always selected (ruling 40). 50 of the 899 ``.py`` files under ``tests/``
  name a process API (`[M]` 2026-11-01, a regex for ``subprocess``,
  ``Popen``, ``multiprocessing``, ``os.fork``, ``posix_spawn`` and
  ``os.system`` with a positive control; the plan's 46 is a narrower regex);
  how many TESTS start one is not counted. The machinery is a recorder of a
  started process's reads and children (#690), the memo's own boundary
  (:ref:`verification-reference-cache-unseen`).

Two further limits are properties of the evidence, not of the store:

* a test reading children from its own temporary memo store (a memo test
  using ``tmp_path``) holds a ``ChildPin`` that ``select`` validates against
  the default store, so its cone is always stale and the test always
  selected (correct, an over-selection);
* the capture's cost on the whole suite has not been measured. The round-2
  extrapolation `[R]` is about 1.4 h of fork overhead on about 20 000 tests,
  paid once per environment and then only for stale tests.


.. _verification-cone-store-gotchas:

Gotchas
=======

* **Run every verb under the suite's interpreter and optimisation level**
  (``python -O``). Both are in the environment key; a ``select`` under plain ``python`` after a
  capture under ``python -O`` lists every test.
* **A capture must not change the system under test.** Anything a worker
  exports into ``os.environ`` reaches every interpreter a test starts (the
  bytecode prefix above flipped a test). Set per-process state inside the
  worker.
* **Code that runs before the recorder is in nobody's cone.** The recorder
  module is loaded by path for that reason; importing anything from
  ``orpheus`` in ``tools/cone_store/__init__.py`` would break it (the
  package imports nothing of ``orpheus`` at import time).
* **A gate that compares two capture routes cannot see a defect in the
  reading both share** (G4 and T1 compare routes through one ``cone_of``):
  the recorder's ground is R1 to R5.
* **Containment is green for every superset.** Never read T1 to T3 as
  evidence that a selection is tight; S1, S7 and P3 bound it on the fixture
  only.


.. _verification-cone-store-history:

History
=======

.. list-table::
   :header-rows: 1
   :widths: 12 58 10 20

   * - When
     - Milestone
     - Issue
     - Where
   * - 2026-11-01
     - **Increment B1**: the capture (grandparent, collecting parent per
       module, child per test, fork guard, fresh fallback), the cone, the
       SQLite store, the two-stage freshness, ``select``, ``status``, the
       environment key; 80 gate ids; the 35-arm battery and the 14-arm
       round after review. Ruling 42: the store predicts full-run
       verdicts, so every conftest is loaded in the grandparent and their
       set, the ``pytest11`` entry points, ``PYTEST_PLUGINS`` and
       ``PYTEST_DISABLE_PLUGIN_AUTOLOAD`` are in the key. The memo gains ``collect_pins`` and the
       ``Refusal`` kinds, keys its code objects by file (ERR-124), builds
       ``trace_call``'s manifest under ``unrecorded()``, leaves
       ``__pycache__`` out of listings, takes a reference moment in
       ``Witness`` and exposes ``fork_blockers``.
     - #358
     - branch ``feature/cone-store``, uncommitted at the time of writing
   * - 2026-11-01
     - **Ruling 38 corrected**: one collecting parent for the suite leaks
       caches and ``sys.modules`` between modules (the fixture's five
       missing nodes); repaired by a grandparent and a parent per module.
     - #358
     - plan ``6bc9a872``
   * - 2026-11-01
     - **Ruling 42**: the store predicts full-run verdicts; a per-module
       verdict is refuted for the question "what does the store promise",
       since no gate reads one.
     - #358
     - plan ``70351e87``
   * - 2026-11-01
     - **Rulings 32 to 41**: capture alone by fork, the environment in the
       key, the store as tooling in SQLite, the row, cone unknown, the
       selection command as the first client.
     - #358
     - plan ``adab473a``, ``49ba6b89``
