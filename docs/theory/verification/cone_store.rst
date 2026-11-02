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
   the mutation battery, the scope boundaries; and, since B2 (#689), the
   import-time statement as the unit of a module's import-time code (the
   two ways a run observes one, the class rules, the projection pin, the
   import-outcome pin, the quiet-section rule, the store's FORMAT). Not
   owned here: the pin kinds, the recorder, validation, the kept reads,
   the fork guard and the refusals (reference_cache.rst, whose sections
   this page cites); the
   dependence DAG of the suite and its rests_on edges
   (principles.rst, vv-test-suite-layout); the V&V markers (harness.rst).
   Plan of record: .claude/plans/verification_tree_refinement.md, from
   "The cone store -- design" to its end (rulings 32 to 44, ruling 38 as
   corrected, ruling 33 as amended by 43). Specification: scratch/verification_tree/cone_store/spec_b1.md;
   reviews qa_b1.md, elegance_b1.md, battery_b1.md, m1_diagnosis.md in the
   same folder; for B2, b2/gate_spec.md, b2/qa_b2_design.md (with the 42
   attack rows in b2/attacks/), b2/qa_b2_build.md and b2/elegance_b2.md.
   Written for increment B1 (landed 7fdc7474), updated for B2 on branch
   feature/cone-store-b2, uncommitted, 2026-11-02, over 549652ed.

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
  (``DefPin``), the import-time statements of a module that the test could
  observe (``StatementsPin``), the outcome of importing the modules it
  imported (``ImportOutcomePin``), a data file's bytes (``DataPin``), what a
  path is (``PresencePin``), a directory's names (``ListingPin``), a memo
  entry read (``ChildPin``) or an installed distribution
  (``DistributionPin``). The interpreter, the platform and the rest of the
  environment are one key per cone, not pins.
* **A module's import-time code is pinned per statement** (B2, #689). A
  run observes a top-level or class-body statement by NAME (code it ran
  reads a name the statement binds) or by EFFECT (the statement does more
  than bind, and runs for every importer). An edit to a constant selects
  the tests whose code reads it, not every importer of its module: on the
  production tree, an edit to one tolerance of ``orpheus/mesh/axis.py``
  selects the 6 tests of ``test_axis_primitive.py`` that ran its reader,
  where module grain selected 9 (`[M]` gate SG4, 2026-11-02;
  :ref:`verification-cone-store-statements`).
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
* **Over-selection is bounded by precision gates on fixtures only.** The
  ground-truth gate is containment (the store's cone holds the alone cone),
  which every superset passes. The precision rows (SG1, SG2 on a fixture;
  SG4 on a clone of the production tree) bound it for one constant each.
  What remains wide is named: the names the grandparent's code read are
  entry names of every test, and every effect statement of an imported
  module is observed by every test that imported it.
* **Evidence.** 92 gate ids in 8 files under ``tests/gates/tools/``
  (89 not ``slow``; `[M]` 2026-11-02, ``pytest --collect-only``), and 91
  memo-level rows in ``tests/gates/numerics/test_traced_memo_statement_pins.py``;
  B1's 35-arm mutation battery with a positive control, every arm reddening
  its target, and a 14-arm round after review, every arm reddening its
  target and no declared-green row; the memo's own suites give the same
  outcomes captured and run plainly (M1, `[M]` 1342 s at B1).


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
   * - ``StatementsPin``
     - a first-party module that ran, by the digest of the ordered list of
       its import-time statements the test could observe, selected on the
       test's entry names (kept once per cone, by their digest, under
       ``.cache/sets``)
     - an observed statement changes, a statement binding an entry name
       is added or removed, the observed statements are reordered, or any
       effect statement changes
   * - ``ImportOutcomePin``
     - the first-party modules the test imported, in import order, with
       the digest of each one's bytes and its dotted name, and the outcome
       of importing them in a fresh interpreter (each module's error and
       warnings)
     - some module's bytes change AND importing the set now gives another
       outcome (an unobserved statement that starts raising or warning)
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
pinned file had at capture (ruling 39). The grain is the memo's: def bodies,
and the import-time statements a test could observe
(:ref:`verification-cone-store-statements`). An edit to a statement no code
of the test reads by name, and that does no more than bind, leaves the
test's cone valid unless it changes what importing the test's modules gives.
The two directions are not symmetric. A
cone larger than the truth costs a test run; a cone smaller than the truth
lets a broken test go unselected. Every rule below is chosen so that a
case it cannot decide makes the cone larger.

A memo entry the test read is a ``ChildPin``, not the defs that generated
it: a warm hit runs no oracle code in the test's process. An edit to the
oracle stales the entry (its own manifest pins the oracle's defs), and the
stale entry stales its reader through the ``ChildPin``. ``select`` judges a
``ChildPin`` by the memo's validation, which imports and runs nothing in its
own process and never regenerates the entry (gate P2 asserts the memo
store's files are byte-identical after ``select``). Validation may start one
process: when a module of an import set has other bytes, the set's import
outcome is measured again in a fresh interpreter
(:ref:`verification-cone-store-import-outcome`).


.. _verification-cone-store-statements:

Import-time statements: what a test observes of a module
========================================================

A first-party module's code is of two kinds. A def body runs only when the
def is called, and the memo pins it by its own digest when it runs
(``DefPin``). Everything else runs when the module is imported: the
module's top-level statements and the statements of every class body. That
code is pinned per statement, by the statements a test could observe, not
per module. The statements, the rule that selects them and the closure are
:mod:`orpheus.numerics._import_time`, a module of the memo's that imports
nothing from ``orpheus`` and is one of the memo's own sources (an edit to it
makes every memo entry stale through ``MemoPin``); the cone store reads the
memo's pins and owns none of this machinery (ruling 15 of the plan).

**A statement** is a top-level statement of a module or, recursively, a
statement of a class body, in source order, a class statement followed by
its body's statements. Its digest is the memo's normalised digest (the
``ast.dump`` without positions, docstrings stripped) of the statement with
every def body read as ``pass``: a def statement is its decorators, its
signature with defaults and annotations, and its name, never its body; a
class statement is its header (decorators, bases, keywords) with its body
read as ``pass``, the body's statements being statements of their own, each
with its class as owner. The index of a source text keeps, per statement,
its digest, the names it binds, the names it reads at import, whether it is
an effect, whether it is an import, whether it is a def, its owner and its
lines (:ref:`verification-reference-cache-source-index`). `[M]` 2026-11-02,
the statement census (``scratch/verification_tree/cone_store/b2/statement_census.py``):
12 339 statements in the 415 files of ``orpheus/**`` and every
``conftest.py`` under ``tests/``; the test modules hold 26 150 more (qa's
census).

A test observes a statement in one of two ways. **By name**: the statement
binds a name (an assignment, a def, a class, an import), and code the test
ran reads that name. **By effect**: the statement does more than bind, so
whatever it does happens for every importer. The pin then digests exactly
the statements observed, and an edit to any other statement leaves the
pin valid; the import outcome (below) is what still catches an unobserved
statement whose edit breaks the import itself.

Observed by name: the entry names, closed over the statements
-------------------------------------------------------------

**The seed.** Every code object the test started contributes its
``co_names`` and the identifier tokens of its string constants, recursing
through tuples and frozensets (``code_names`` and ``string_names`` in
:mod:`orpheus.numerics._import_time`):

* ``co_names`` holds every name a code object reads as a global or an
  attribute (``LOAD_GLOBAL``, ``LOAD_ATTR``) and every global it stores
  (``STORE_GLOBAL``, so a ``global X; X = …`` inside a def reads as the
  name ``X``);
* the identifier tokens of string constants catch a name read through a
  string: ``getattr(x, "name")``, ``hasattr``, ``attrgetter("inner.val")``
  (both tokens), a format string ``"{0.label}"``, a ``match`` class
  pattern's keywords, iteration over a set of names (a set literal is a
  ``frozenset`` constant, hence the recursion through frozensets);
* standard-library, third-party and generated code contribute too: a
  dataclass's generated ``__init__`` names every field.

Two kinds of code add more, and one kind adds nothing:

* **A first-party def that ran adds the parts of its qualified name**
  (``Class.method`` gives ``Class`` and ``method``; ``<locals>`` is
  dropped). Pytest finds a test function by iterating its module's
  namespace, so a test's own name is in no code object's ``co_names``, and
  without this rule the test's def statement, whose decorators hold its
  ``parametrize`` list, would never be selected (qa's row S2a: 371
  module-level names in 183 of 769 test modules are read only by a test or
  fixture decorator, `[M]` 2026-11-02; and S2b″, a ``fixture(params=…)``
  constant, the fixture being found by the test's argument name, which is
  in ``co_varnames``, not ``co_names``).
* **A first-party lambda or generator at module or class level that ran
  after import** adds its own names, and the names bound by the statement
  whose lines hold it (``holding``), so the statement that defines it is
  selected.
* **First-party module-level code adds nothing**: a module's ``<module>``
  code, a class body, a signature scope (``__annotate__``, a PEP 695
  scope). Its statements are read statically, one by one. Counted as
  executed code, a module's ``<module>`` reads nearly every name its module
  binds, so every statement of every imported module would be observed and
  the statement grain would gain nothing over the module grain; dropped
  without a replacement, every read made at import (a default argument, a
  derived constant ``B = A * 2``, a re-export) would be lost. The static
  reads below are the replacement.

**What a statement reads and binds.** A statement's reads are the ``Name``
loads, the ``Attribute`` names and the string-constant tokens of what it
evaluates at import, with a def's decorators, defaults and annotations
(annotations are lazy under PEP 649, but reading them evaluates them, so
they are reads), a class's decorators, bases and keywords, and lambda
bodies; an ``ImportFrom`` reads the names it takes from its module. Its
binds are the names it stores or deletes in the namespace that holds it (a
``del`` target counts as a binding), a def's or class's name, an import's
alias (the first component of ``import a.b``), an ``except … as`` name, a
``match`` capture, and a walrus target, which binds in the enclosing
namespace even inside a comprehension; a name stored inside a lambda or a
comprehension binds nothing at module or class level.

**The closure.** The entry names E′ are the least fixed point of: the
seed; plus the reads of every statement selected under the current names,
across every first-party module the test ran (``closure``). An import
adds its reads only when it is selected by name, that is, when the test
reads a name the import binds: a re-export is then followed to the module
that defines the name. An import selected only because it is an effect
(below) adds nothing to the closure, or every name every import mentions
would be read. Names are not qualified by module, so every module binding
a name the test read is selected: a collision over-selects, and never
under-selects.

Observed by effect
------------------

A statement is an **effect**, and every test that imported its module
observes it, unless it is a pure binding. The predicate, as
``_import_time._effect`` spells it, over the nodes the statement evaluates
at import (never a def's or a lambda's body, never an annotation):

* an import (the imported module's own effects are its own pins);
* a bare expression that is not a constant (a call that binds nothing);
* a ``del``, ``raise``, ``assert`` or ``global`` anywhere in the statement,
  a nested import, and a store into a subscript or an attribute;
* a call whose callee is not pure, and a call to a pure higher-order callee
  (``map``, ``filter``, ``sorted``, ``min``, ``max``) handed a function
  whose root is not trusted (``list(map(register, PLUGINS))``; a lambda
  counts as untrusted);
* a def or class with a decorator that is not a pure callee;
* a class with a base or a keyword whose root is not trusted: creating it
  can run a first-party ``__init_subclass__``, a ``__set_name__`` or a
  metaclass.

A callee is **pure** when its last name is on a short list of callees that
only return a value (builtins such as ``frozenset`` and ``sorted``,
``dataclass`` and ``field``, ``numpy`` constructors such as ``zeros`` and
``linspace``, the SymPy and mpmath number constructors, ``pytest``'s
``mark``, ``parametrize``, ``fixture`` and ``param``) and its root is
trusted. A root is trusted when it is a builtin the module never rebinds,
or a name bound only by imports from a short list of modules
(``builtins``, ``dataclasses``, ``functools``, ``typing``, ``collections``,
``enum``, ``abc``, ``pathlib``, ``re``, ``math``, ``operator``, ``types``,
``fractions``, ``decimal``, ``numpy``, ``pytest``, ``sympy``, ``mpmath``,
and ``typing_extensions``). So a first-party ``cache`` or ``field`` is an
effect, a bare ``class K(Enum)`` or a class with metaclass ``ABCMeta``
imported from ``abc`` is not, and a first-party mixin beside ``Enum`` or
``Protocol`` is (`[M]` qa's build review, one ``_effect`` call each).
``type`` is not on the list: ``type("Bar", (Base,), {})`` creates a class,
and a first-party ``Base`` with a registering ``__init_subclass__`` is an
effect like ``class Bar(Base)`` (qa's row T1).

**Why the complement and not "calls first-party code".** The B2 design
selected as effects only the statements whose import-time evaluation calls
first-party code. qa's adversarial pass found seven shapes it misses, each
realised on a fixture: a third-party effect (``np.seterr``,
``warnings.filterwarnings``), a builtin method mutating a module object
(``DEFAULTS.setdefault``), a first-party instance bound by assignment
(``registry.add``), a store or call nested in a compound statement (a
module-level ``for`` filling a table), class creation, which runs a
first-party ``__init_subclass__`` with no ``Call`` node, and an import that
registers. So an effect is anything that is not a pure binding. The price
`[M]` (qa's model, ``attacks/repair_cost.py``): 2 119 of the census's
12 339 statements are effects, 1 231 of them first-party imports, against
248 under the design's predicate; the other 10 220 are selected by name.

The class rules
---------------

Within a class, selection follows four rules, applied with the two above
to a fixed point (``selection``):

* a class header is selected with any statement of its body (a method
  observed by name selects its class);
* a selected header selects every non-def statement of its body: a field
  iterated by ``dataclasses.fields``, an enum member iterated with its
  class, an attribute a first-party ``__init_subclass__`` reads, none of
  them read by name;
* a selected header selects every statement of its body binding a dunder,
  a def included (``__hash__ = None``, ``__match_args__`` and
  ``__array_priority__`` are read by the interpreter or by numpy, never by
  Python code by name; a dunder method's default is evaluated at import);
* a class statement is itself an effect when a base or a keyword is not
  trusted (above), and then its header is observed by every importer.

Their populations `[M]` (qa's census of ``orpheus/``): 2 802 of the 12 339
statements are non-def class-body statements; 462 dataclasses, with 104
field-iterating readers, among them the content identity
(``orpheus/values/content.py``); 22 enums; 159 classes with a first-party
base in another module.

The pin: one projection per module, on one kept name set
--------------------------------------------------------

``StatementsPin(path, names, digest)``: ``names`` is the digest of the
test's whole E′, kept once under ``.cache/sets/<digest>.json`` (canonical
JSON, named by the digest of its bytes); ``digest`` is the module's
**projection** on E′, the digest of the ordered digests of the module's
selected statements (``_SourceIndex.projection``). Every module whose
code the test ran gets one, and all of a test's ``StatementsPin`` rows name
the same set. Validation reads the kept set and re-runs the selection on
the module's CURRENT statements with the STORED names, so the pin is
stale when:

* a statement the test observes changes;
* a statement binding one of the stored names is added (a shadowing def, a
  new dataclass field of a selected class) or removed;
* the observed statements are reordered;
* any effect statement changes, or one is added or removed.

An edit to a statement the test does not observe leaves the pin valid.
`[M]` qa's build review (``build/rows_reproj.py``): appending ``VALUE = 9``,
a new statement binding a stored name, makes the pin stale; appending an
unread statement leaves it valid and the result unchanged.

**Why the whole E′, and not the names occurring in each module.** The
design stored, per module, E′ restricted to the identifiers of that
module's source at capture. A name new to a module is then never in its
stored set, and four realised shapes read a name a module did not hold at
capture: a subclass in another module gaining an override (``tol = 1e-3``
on ``Child(Base)``); a package whose ``__getattr__`` served ``FOO`` gaining
a real ``FOO``; ``from pkg import thing`` (a submodule) shadowed by a new
``class thing`` in the package's ``__init__``; and an unguarded attribute
read that raised at capture and finds the name after the edit. Storing E′
whole once per test (ruling 43, R4) closes the four (qa's rows S3a, S3d,
S3e, S3f), and costs one set per test rather than one per module.

**Why one pin per module rather than one per statement.** The store judges
content pins file by file, and a file whose identity is unchanged skips
every pin on it (:ref:`verification-cone-store-freshness`). One pin per
module keeps that stage at file grain and keeps the cone's rows per file;
the statement grain lives inside the projection.

**The kept set is content-addressed.** Its name is the digest of its bytes,
so any process may write it and no process can write a wrong one under a
name: two writers write the same bytes, atomically. A set file whose bytes
do not digest to its name is written again by the next generation that
keeps it, and a read that finds it missing or corrupt is not remembered for
the process, so a later generation can repair it. A pin whose set is
missing is stale ("the entry names … are not kept"), which over-selects.

**The size, on one real run** `[M]` (qa's build review,
``build/warm_count.py``: a trace of ``import orpheus.homogeneous``): 10 580
entry names; 85 modules in the import set; 93 ``StatementsPin``;
1 794 of the 2 348 import-time statements of those modules selected
(76 %). The selection decomposes as: 869 effects; 1 035 with the class
rules applied to the effects; 1 411 with the closure of an empty seed;
383 more from the names of the third-party and frozen code the run
started. That is the statement fraction of one import-only run, not the
fraction of memo entries an edit leaves warm (below).

.. _verification-cone-store-import-outcome:

The import outcome: an unobserved statement that breaks the import
------------------------------------------------------------------

A projection is blind to one class of edit by construction: an edit to a
statement the test does not observe that makes the import itself raise or
warn. Deleting a def that no test calls while an unread table still names
it (``_LEGACY = {"old": old}``) makes the module's import raise
``NameError``, and every test that imports the module changes verdict
while every projection is valid (qa's row S2k); the same across modules,
``from .b import old`` in an unread binding of ``a.py`` when ``b.py``
deletes ``old`` (S2k′); and an unobserved statement that starts warning
under a session ``error`` filter (S2l, the shape of the suite's
``DriftWarning`` filter). Every unselected statement is exposed.

``ImportOutcomePin(modules, outcome)`` closes it (ruling 43, R3 amended).
At capture, the import set is every first-party file whose ``<module>``
code started, in the order it started, each as (pin path, digest of the
bytes it was compiled from, the dotted name it was imported by,
``__spec__.name``), kept under ``.cache/sets`` by the digest of the list.
The outcome is the digest of what a fresh interpreter reports when it
imports every module of the set by its dotted name, in that order: per
module, the error that stopped its import (or none) and the sorted set of
warnings it emitted, with every warning shown (``simplefilter("always")``).
Validation compares each row's bytes digest with the file's current bytes;
when all are equal the pin is valid with no measurement, and otherwise the
set is measured again as it now is and the two outcomes are compared.

**The instrument.** ``_import_outcome`` in
:mod:`orpheus.numerics.traced_memo` runs ``python -B -c`` with a program
held as a string in the memo's source, the set's dotted names and their
roots on its standard input, under the generation's environment and inside
the memo's ``spawning()``, so the measurement is never recorded as the
test's own process start. ``-B`` keeps it from writing bytecode beside the
sources. A module's root is the directory that, joined with the path of
its dotted name, gives the real file the run executed (``src`` for
``src/libmod.py`` imported as ``libmod``; ``_module_roots``), derived from
that real file and never from the row's pin path, which is relative to the
first ``sys.path`` entry holding the file. The measuring interpreter puts
the roots ahead of its own default path, ordered as the importing
process's ``sys.path`` orders them: two roots can hold one top-level name
(the repository's own ``tests`` package and a test tree's), and in the
set's import order the repository's came first and shadowed the tree's,
so every test module of the tree became unmeasurable (the test-architect's
SG0, 2026-11-02, after the first form of this repair). The measuring process's
``sys.path`` takes no part: it is the run's at capture and the validator's
afterwards, and a module the run found through pytest's ``pythonpath``
setting was not importable by ``select``, so a comment edit to it gave
another outcome and selected every test that imported it (`[M]` the
test-architect's gate SG7, 2026-11-02; its import-raises leg was stale for
``ModuleNotFoundError``, the wrong reason). The result is kept under
``.cache/import_outcomes/<key>.json``, where the key is the digest of the
memo's own source digest (the one a source index is named by, so an edit to
the measuring program renames every outcome), the interpreter, the
environment, the roots and the rows: one checkout state is measured once,
and every entry or cone sharing a set shares its measurement. A file whose bytes changed between the run and its manifest
records the outcome ``_UNMEASURED``, which no measurement equals: the pin is
stale until the file holds the bytes that ran again.

Four alternatives were refused, each with its reason:

* ``-W error`` (qa's proposal R3W): a warning already emitted at
  capture would make the check fail on every later validation;
* re-importing only the modules whose bytes changed: the module that
  breaks need not be the one that changed (S2k′: ``a.py`` raises because
  ``b.py`` changed; `[M]` qa's model, the row lands under that rule);
* ``-m`` on a module of the memo's: running a module of
  ``orpheus.numerics`` imports ``orpheus/numerics/__init__.py``, which
  imports most of ``numerics/`` into the measuring process and changes what
  it measures (ruling 44);
* loading a module from its file path: outside its package a module's
  relative imports fail, and a dataclass under ``from __future__ import
  annotations`` fails looking its module up in ``sys.modules``, both
  before the body runs, so the outcome could not move. `[M]` qa's build
  review (finding F1, ``build/census_standalone.py``): by absolute path,
  256 of the 404 modules under ``orpheus/`` raise before their body runs
  (131 dataclass, 118 relative import, 7 other); by dotted name, 0 of 404.
  The route was live for any generation started outside pytest and the
  docs build.

F1 is closed in two parts, both refusals rather than measurements that
cannot move:

* the measurement imports each module by the dotted name recorded at
  capture (``__spec__.name``), and a first-party module the run imported
  that no ``sys.modules`` entry names is refused (``UnnamedModule``, which
  makes a cone unknown);
* a module the run imported that a fresh interpreter cannot import (the
  run reached it through import machinery of its own process, a finder it
  installed in-process, say) would fail at capture and at validation
  alike, so its outcome could not move either. Every module of the set
  imported in the run, so any error in the capture's measurement means the
  instrument does not reach that module as the run did, and capture
  refuses the pin (``UnmeasurableImport(name, error)``, from
  ``_measured_import_outcome``; ``_import_outcome`` returns the digest with
  each module's error). An editable install's finder is not such a case: it
  is installed from a ``.pth`` file at interpreter start, so the measuring
  interpreter has it too.

**The cost, and #691.** A pin with no file is judged on every ``select``.
With every module's bytes unchanged that is one digest per module path.
After an edit to a broadly imported module, every distinct set is measured
again once: `[M]` qa's build review (``build/store_sets.py``) on the 73
tests of ``tests/gates/values``, 16 distinct sets of 249 to 252 modules
each, every one containing ``orpheus/numerics/space.py``, took 12.2 s to
re-measure after an edit to that file (0.76 s a set; one cold measurement
of an import-only run's set, 0.65 s). `[R]` Over the 20 661 collected
tests that is about 4 500 sets and about 57 minutes for the first
``select`` after such an edit. And the set keeps the absolute paths of
modules a test writes and imports under pytest's ``tmp_path`` (`[M]` 14 of
the 73 tests): a renamed or rotated base temporary directory reads those
modules as absent, and the test is selected. Both are #691 (measure the
outcome per module once over the union of the stored sets; filter the
rows to first-party files, relative to the repository), whose candidate
repair owes a soundness argument first, because a module's import would
then be measured in the union's order, not in the run's.

Code the memo runs for itself is no run's
-----------------------------------------

The recorder (``Recording``, a ``sys.monitoring`` ``PY_START`` callback)
records each code object the first time it starts and then disables its
event. Code that starts inside the memo's ``unrecorded()`` section (the
memo's own reads, its pin collection, the fork template's own work) is
**not recorded and not disabled**: its event stays enabled, so the code is
recorded when the run itself first starts it. And the fork template runs
its own pre-fork work (the snapshot of the recording, each job's checks, the
fork and the child's restore) inside ``unrecorded()``, so a forked child
records what a fresh process records and nothing of how it was born.

The rule matters because every started code object now contributes entry
names. Before it, the template's own work was recorded in a forked child
and not in a fresh process. `[M]` the test-architect, 2026-11-02, the two
kept name sets of gate E1's fixture
(``test_traced_memo_template.py::test_e1_a_forked_childs_manifest_is_a_fresh_processs``):
22 073 names in the forked child against 22 054 in the fresh process; the
19 extra names (``_take``, ``snapshot``, ``getpid``, ``makedirs``,
``_copy_atomic_types``, …) came from ``Recording.restore`` and the copy
and directory calls the memo makes in the child. Every ``StatementsPin``
of an entry names the one set, so all of them differed: 671 pin rows.
Gate E8 (a nested miss against a direct miss) diverged the same way. At B1
the module grain did not depend on the names, so both gates passed. Their
run after the repair is pending at the time of writing.

.. dropdown:: First got wrong: "a name collision can only over-select"
   :color: muted

   Ruling 33 chose names as the unit and argued that, since names are not
   qualified by module, a collision can only over-select. The B2 design
   built on it: the entry names from executed code, module-level code read
   per statement, an effect being a statement whose import-time evaluation
   calls first-party code, and each module's pin projected on the names
   occurring in its own source. qa modelled the design as written and
   attacked it with 42 fixtures (``attacks/attacks.py``, two of them
   positive controls): `[M]` 24 of 42 rows under-selected, every pin valid
   while the run's result changed. The sub-claim is [REFUTED 2026-11-02]
   for the question "is the name rule sound"; the fact it establishes is
   that an import-time statement can change a run without binding any name
   the run reads. The shapes, with their populations in the tree: a
   constant read only by a ``parametrize`` or ``fixture(params=…)``
   decorator (371 names in 183 test modules); an unselected statement that
   starts raising at import (every unselected statement); a subclass
   override, a new dataclass field or enum member read by iteration (159
   cross-module subclasses, 462 dataclasses, 22 enums); a subclass added to
   an ``__init_subclass__`` registry
   (``orpheus/sn/operators/loss_kernel_gauge.py`` iterates one); a table
   read only by a module-level loop (``orpheus/data/group_structures/wims.py``,
   saved on the real tree only by a collision on the loop variables);
   dotted and set-literal strings. Ruling 43 repaired it under ruling 33's
   intent with six amendments, each closing named rows: R1 (a def that ran
   selects its qualified name's parts), R2D (an effect is anything not a
   pure binding), R3 (the import outcome), R4 (E′ stored whole), R5 (a
   selected header selects its body's non-def statements) and R6 (string
   tokens recursing through frozensets). Together they leave 0 of 42 rows
   landing in qa's model and keep control 2's benefit; dropping any one
   reddens rows (R1: F7; R2D: 9 rows; R3: 4; R4: 4; R5: 6; R6: N4 and N5;
   `[M]` the test-architect's ablation on qa's model, ``ta/ablate.log``).
   Refuted attacks, each for a structural reason: a walrus, a star target,
   a PEP 695 ``type`` and an augmented assignment are ``Store`` targets;
   ``global`` in a def is in ``co_names``; ``__slots__`` is a dunder;
   a module-level lambda is selected by line; a re-export and
   ``import a.b as c`` follow unqualified names; a field reorder moves the
   ordered digest.

.. dropdown:: First got wrong: "the ModulePin class stays, so old rows still parse"
   :color: muted

   Step 3 of the B2 design kept ``ModulePin`` for one cycle so that rows
   written before B2 would still decode. The build deleted the class, and
   the elegance review found that a store written by B1 then crashed
   ``select`` (``KeyError: 'ModulePin'``). Ruling 44 retired the module
   skeleton and ``ModulePin`` with no shim: every stored memo entry is stale
   anyway, because the memo's own source changed, and an old memo row reads
   as ``Corrupt`` and regenerates; the cone store's rows carry
   ``FORMAT = 2`` and a store of another form is emptied when opened. An
   always-stale ``RetiredPin`` was refused: it would keep a dead kind in
   the ``Pin`` union, which the format version makes unnecessary.

The evidence for B2
-------------------

* **Memo level**,
  ``tests/gates/numerics/test_traced_memo_statement_pins.py`` (91 rows,
  `[M]` ``pytest --collect-only``, 2026-11-02), one small package per case
  (``_statement_pin_cases.py``; qa's 42 rows transcribed by script into
  ``_statement_pin_qa_rows.py``, never by hand). Each row traces a call in
  one fresh process, edits, and validates in a second. A precision row
  asserts that ``validate`` returns no reason at all after an unobserved
  edit (an unread constant, an unread derived constant, a reorder of
  unselected statements, an unread def's default); a recall row asserts
  the manifest goes stale through a statement or outcome pin, that no
  ``DefPin`` saw the edit, and that the edit changes the run's answer (X1);
  an under row pins a declared boundary (valid while observed); an over
  row a declared over-selection. The positive control on B1's code: with
  ``ModulePin.stale`` disabled in the validating and selecting processes,
  81 of the 92 rows of the two files then present reddened, and the 11 that
  did not are exactly the six precision rows, the two boundary rows, SG0
  and the two population rows (`[M]` ``ta/baseline_arm_blind_module.log``).
* **The memo's own witness inverted.** Gate M2.1's row ``constant-unread``
  (an edit to a module constant no traced code reads) asserted a stale
  entry at B1, a declared spurious miss of the module skeleton; at B2 it
  asserts the entry stays valid.
* **The store**, ``tests/gates/tools/test_cone_store_statement_grain.py``:
  SG0 (no edit selects nothing; pytest collects exactly the hand-written
  ids); SG1 and SG2, control 1 exact on a fixture (red at B1: the sibling
  and the neutral test selected too); SG3, recall at the store (an effect
  statement, a renamed re-export, a default argument, qa's S2a and S2b″);
  SG4 (``slow``), control 1 exact on the production tree; SG5, a store in
  B1's form. **SG4**: a copy-on-write clone of ``orpheus/``, ``tests/`` and
  ``tools/``, the tests of ``tests/gates/sn/primitives/test_axis_primitive.py``
  captured, ``_OUTFLOW_COSINE_TOL`` of ``orpheus/mesh/axis.py`` edited in the
  clone, and ``select`` required to equal the readers found by an
  independent monitor (each test alone in a fresh interpreter, a
  ``sys.monitoring`` watch on ``face_outflow_ordinates``). `[M]` 2026-11-02:
  6 tests selected, the 6 readers; B1 selected 9 of 9 (44 s, red).
* **The memo's warm count** (gate 4 of the plan, the precision B2 exists
  for). Before, `[M]` the test-architect on ``e9df7ed9`` with
  ``ta/warm_count.py`` over ``orpheus/numerics/metric.py``'s
  ``_DENSE_METRIC_PENROSE_RTOL``: 7 835 of the 7 835 readable entries of the
  store go stale, since every entry pins ``metric.py`` by module (33
  modules are pinned by every entry; 691 more entries are unreadable rows
  written before ``98dc5242`` and regenerate as corrupt). After: not yet
  measurable, because no stored entry carries a ``StatementsPin`` until it
  regenerates under B2. `[R]` 0 of the 7 835 entries entered
  ``DenseMetric.__post_init__``, the only def that names the constant, so
  the after count is expected near 0, less the name collisions.


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
  importer, so its module code is attributed to whichever test module's
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
every test is refuted for the same question: any edit to the import-time
code of any module some collection ran would select the whole suite.

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

**Why the two routes record the same names.** At B2 a test's cone depends
on the names its run reads (:ref:`verification-cone-store-statements`),
so any code the capture
itself runs on one route and not the other changes the test's rows. The
test-architect's route probe (G4: one test captured by fork and by
``--route fresh``, its entry names compared) found five such sources in
turn, and each is now run unrecorded or not run at all:

* each role's bookkeeping (the fork guard, the forks, the after-fork
  hooks), run under ``unrecorded()``;
* the terminal reporter, whose progress and width code ran inside each
  test's protocol differently per route (``-p no:terminal``; a capture's
  verdict comes from the reports, never from its output);
* bytecode writes, which a fork child made for modules its parent had
  already compiled (``sys.dont_write_bytecode``);
* the memo's query of numpy's BLAS (``_blas``), whose numpy call is cached
  per process: the fork route's guard made it before any test, the fresh
  route inside the test. The cached call is numpy's own
  ``show_config`` (``_NUMPY_BUILD_CONFIG``, a ``functools.partial``, so no
  first-party code is cached) and runs unrecorded; ``_blas`` itself reads
  the name out of it on every call, recorded, so an edit to either is
  pinned (qa's final pass, F3: a cached first-party body was pinned by no
  cone);
* pytest's sort of a module's items by fixture scope (``reorder_items``),
  which compares items only when there are several: a collecting parent
  sorts the whole module, a fresh child its one test. A conftest's own
  ``pytest_collection_modifyitems`` stays recorded.

A code object started while the memo's ``unrecorded()`` is on is skipped
by the recorder, not disabled, so the same code run later by the test is
still recorded.

**The capture's own code is in no cone.** The worker module runs in every
capture and in no full run, so it is no test's dependency. `[M]` 2026-11-02
(a capture of ``tests/gates/test_layer_imports.py``, 533 tests): six of its
defs (the hooks pytest calls inside a recorded stretch, and the wrappers
that switch the recording) were pinned in every cone, and so were the
import-time statements of ``cone.py``, ``pins.py``, ``store.py`` and
``worker.py``, so an edit to the store's own code selected the whole suite
(qa's and the elegance-enforcer's final pass on B2, F5). The worker removes
its own module's code objects from the recording before it reads it, and
imports the modules that read the cone only after the recording stopped. A
test of the cone store imports those modules itself, during its own run,
where they are recorded as any other code. Turning the recording off and on
is one switch, the boot module's ``unrecorded()`` and its twin
``recorded()``, which every wrapper of the worker goes through.

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

* ``DefPin`` and ``StatementsPin`` of first-party code: any source under
  the repository that is not the interpreter's or an installed
  distribution's (a virtual environment inside the checkout) and not the
  bytecode cache. A test's own body and its helpers under ``tests/`` are
  first-party, so an edit to a test selects it at def grain, and so is
  ``tools/``. B1 kept only ``<repo>/<package>/`` and ``<repo>/tests/``, so
  an edit to ``tools/cone_store/`` selected none of the tests under
  ``tests/gates/tools`` that ran it, an under-selection found in B2 when an
  import-outcome row named a ``tools/`` module the store held no pin on;
* the ``ImportOutcomePin`` as it comes (its rows are not filtered to
  first-party files: :ref:`verification-cone-store-import-outcome`). At
  ``select``, stage one judges it by the rows whose file the store holds a
  pin on; a row it holds none on (a module a test wrote under
  ``tmp_path``, code of another checkout on ``sys.path``) is read only when
  the set is measured again (#691);
* ``DataPin``, ``PresencePin`` and ``ListingPin`` of paths under the
  repository, except under the bytecode cache directory and the directories
  above it that the capture itself creates;
* ``ChildPin`` and ``DistributionPin`` as they come.

It drops, by decision, the interpreter, environment, working-directory and
memo-source pins: the environment key covers the first two, and the memo's
own sources are first-party files whose defs and import-time statements are
already pinned. The paths of every pin that names a path are stored relative
to the repository, so such pins survive a move of the checkout (the rows
inside an ``ImportOutcomePin``'s kept set are the exception: they keep the
memo's pin paths, absolute for a module found outside ``sys.path``, and a
moved checkout re-measures the outcome, an over-selection);
``tools/cone_store/pins.py`` is the one
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
   * - ``NoRecordedBytes``, ``SourceGone``, ``NoDefAtLine``, ``UnnamedModule``
     - when the file is first-party; dropped otherwise (``UnnamedModule``:
       the test imported a first-party module no ``sys.modules`` entry
       names, so its import outcome cannot be measured again)
   * - ``AmbiguousOwner``, ``UnownedSitePackage``, ``UninstalledDistribution``,
       ``UnmeasurableImport``
     - always (``UnmeasurableImport``: a module the test imported that a
       fresh interpreter does not import, so its import outcome cannot be
       measured again)
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
       when it is a content pin (``DefPin``, ``StatementsPin``, ``DataPin``)
   * - ``files``
     - one row per judged file: its identity at capture and whether that
       identity is trusted
   * - ``cones``
     - (test id, pin id): the cone of each test
   * - ``meta``
     - the store's ``format``, the grandparent's pin count and the number
       of tests the last capture pruned
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

**The store's form is versioned.** ``tools/cone_store/store.py`` declares
``FORMAT = 2`` (since B2) and keeps it in ``meta``. A store opened with
another ``format``, or with none (every store B1 wrote), is emptied in one
transaction when it is opened, and every test is then new: selected, and
captured again. B1's rows hold ``ModulePin``, a kind the memo no longer
decodes (`[M]` the elegance review: ``Manifest.from_json`` of a
``["ModulePin", …]`` row raises ``KeyError: 'ModulePin'``, which ``select``
did not catch). Emptying is safe in the direction that matters: it selects
everything. Gate SG5 writes a store in B1's form and requires ``select`` to
list every collected test; dropping the format check reddens it with that
``KeyError``.


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
``DistributionPin``, ``ImportOutcomePin``) are judged directly, on every
``select``. Each row of an ``ImportOutcomePin``'s set is judged by a
digest: for a repository file whose identity passes stage one, the digest
the store keeps beside that identity (the digest of the bytes the identity
was taken on, at capture or at the last clean ``select``), so an unchanged
set reads nothing; a file that fails stage one is read and digested (F1
and F4 found every module of the set read for one edited file). A row that
is not the repository's (an absolute path: a module a test wrote under
``tmp_path``, code of another checkout) keeps its own digest and is
measured only when the set is (#691). When any digest differs from the
row's, the set's outcome is measured again
(:ref:`verification-cone-store-import-outcome`). At capture the set's rows
are re-kept in the store's spelling (``relocated_rows``): a first-party
module's row names the same repository-relative file as its
``StatementsPin``, so the stage-one identity it is judged by is that file's.

.. dropdown:: First got wrong: "a file that passes stage one holds the bytes the import set was measured on"
   :color: muted

   The first form vouched a row's kept (capture) digest whenever its file
   passed stage one. But a file's identity is refreshed after a ``select``
   that finds every pin ON the file clean, and the import set names the
   file from outside it. qa's final pass on B2 (F1): delete ``def g`` while
   ``_REG = [g]`` stays at import time. No def pin moves, the import now
   raises ``NameError``, and the first ``select`` finds the import set
   stale; it also refreshes the file's identity, so the second ``select``
   vouched the capture's digest for the edited bytes and read the set as
   fresh (`[M]` qa's probe ``p1_refresh_evaporates.py``: ``{}`` on the
   second ``select``). An identity vouches for the bytes it was taken on,
   never for the capture's, so the store keeps that digest with it.

The memo resolves a ``DefPin`` or ``StatementsPin`` path through ``sys.path``,
so ``select`` puts the repository first on ``sys.path`` while it judges.
A ``StatementsPin`` is a content pin: stage one skips it with its file, and
stage two re-projects the file's current statements on the pin's kept
entry names.

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
   gives (the def or module gone, the digest moved, an observed import-time
   statement changed, the import set's outcome changed, the file's bytes
   changed, the child no longer valid, the path's presence changed).

A test that is no longer collected is not listed. A test parametrised anew
is new under its new ids. A decorator is part of its def's statement, so a
decorator edit selects the tests that observe that statement, and every
importer when the decorator is not a pure callee (an effect). With no edit
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
       pinned where it runs, by def and by import-time statement. Recorded, the whole-file
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

1. **A constant edit selects its readers** (P1, with C3, and exactly SG1,
   SG2 and SG4): the edit to ``TOL`` in ``pkg/consts.py`` selects the tests
   whose code reads it. On P1's tree every importer of the constant's
   module reads it, so P1 and C3 cannot tell statement grain from module
   grain; SG1 and SG2 can, on a fixture where three tests import one module
   and each reads a different constant or none (red at B1: all three
   selected), and SG4 on a clone of the production tree
   (:ref:`verification-cone-store-statements`).
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
``foundation`` (`[M]` 2026-11-02, ``pytest --collect-only``: 92 ids, 89 not
``slow``, every one ``foundation``):

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
   * - ``test_cone_store_statement_grain.py``
     - SG0 to SG5
     - B2's statement grain at the store: no edit selects nothing (SG0);
       control 1 exact on a fixture (SG1, SG2) and on a clone of the
       production tree (SG4, ``slow``); recall at the store for an effect,
       a re-export, a default argument and qa's S2a and S2b″ (SG3); a store
       in B1's form is emptied and every test is new (SG5)
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
(``slow``): 1 passed, 1342 s. The directory collected 181 non-slow ids
with S18 to S22 at B1.
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
  loader (#688), a scope boundary ruled by the orchestrator on 2026-11-01.
* **A name the statement rule cannot see is not seen.** The import-time
  statements a test observes are found by the names its executed code
  could read (``co_names`` and the identifier tokens of its strings), closed
  over the statements (:ref:`verification-cone-store-statements`). The
  machinery that would lift this is a recorder of name lookups in module
  and class namespaces (ruling 43, after qa's 42 attack rows). Not seen,
  each with its population:

  - a name built from data and absent from every string constant of the
    code that ran (gate row U1, an under row: a computed ``getattr`` name
    on a module). The census `[M]` 2026-11-02 counts 36 computed-name
    ``getattr``, 2 ``hasattr`` and 3 ``setattr`` sites and 2
    ``import_module`` calls in ``orpheus/`` and the conftests, and 235
    computed-name ``getattr``, ``hasattr`` and ``setattr`` sites in
    ``tests/``, whose reader is often the test; the sites whose names come
    from a constant tuple are covered by the string tokens, the rest were
    not read one by one;
  - iteration over a module namespace, ``vars(module)`` or a module's
    ``__dict__``, meeting a name new to the module (gate row U3, an under
    row; 4 ``vars()`` and 1 ``.__dict__`` in ``orpheus/`` and the
    conftests);
  - C code reading a non-dunder class attribute of a class whose bases
    and keywords are all trusted roots; a class with an untrusted base,
    a ``ctypes.Structure`` subclass among them, is an effect, and its
    ``_fields_`` is observed (gate row U4, a recall row).

  The guarded probe of a name new to a module (``getattr(module, name,
  default)``), a boundary of the B2 design, is closed by storing the entry
  names whole (gate row U2, now a recall row).
* **Over-selection is bounded but not small.** The names the grandparent's
  code read (pytest's, the conftests' and the 91 ``orpheus`` modules the
  root conftest imports, `[M]` 2026-11-01) are entry names of every test:
  `[M]` qa's design pass measured 24 411 entry names for one test under a
  pytest run against 245 for the same code run plainly. Every effect
  statement of every module a test imported is observed by it, so an edit
  to an import, a registering decorator or a module-level call in one of
  those 91 modules still selects the whole suite. An edit to any byte of a
  module in a test's import set re-measures that set's import outcome, and
  a set holding a ``tmp_path`` module is selected by a rotation of pytest's
  base temporary directory (#691). A new file selects every test that
  probed or listed its directory.
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
  evidence that a selection is tight; S1, S7, P3, SG1 and SG2 bound it on
  fixtures, and SG4 for one constant of the production tree.
* **The import outcome is measured without** ``-O``, while the suite runs
  under ``python -O``. A module whose import-time behaviour differs under
  ``-O`` (an ``if __debug__:`` arm, a module-level ``assert``) is measured
  as it imports without it. `[M]` qa's build review: by dotted name, 0 of
  the 404 ``orpheus`` modules raise on import without ``-O``, so no
  module-level ``assert`` blinds it today.
* **A cone's statement pins depend on a directory the memo owns.** Every
  ``StatementsPin`` and ``ImportOutcomePin`` names a set kept under the
  memo's ``.cache/sets``. Removing that directory makes every such pin
  stale (the set is "not kept"), so the next ``select`` lists every test:
  safe, and a full re-capture's worth of cost.
* **An unread constant is no longer a spurious miss.** A gate or a
  docstring that expects an edit to an unread module constant to make an
  entry stale, or to select every importer, describes the module skeleton,
  retired at B2; such rows were re-posed (M2.1 ``constant-unread``, C4 and
  C5 of the collection gates), and a new one fails.


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
   * - 2026-11-02
     - **Increment B2: import-time statements as the unit** (ruling 33 as
       amended by ruling 43). The module skeleton and ``ModulePin`` are
       retired with no shim (ruling 44); ``StatementsPin`` (one projection
       per module on the run's kept entry names), ``ImportOutcomePin`` (the
       import set re-imported by dotted name in a fresh interpreter when a
       byte of it changes) and the refusals ``UnnamedModule`` and ``UnmeasurableImport`` replace them;
       the statement model is ``orpheus/numerics/_import_time.py``; the
       recorder ignores code started under ``unrecorded()`` and the fork
       template's own work runs there; the store's rows carry
       ``FORMAT = 2``. Gates SG0 to SG5 and the memo's 91 statement-pin
       rows; control 1 exact on the production tree (SG4: 6 selected where
       B1 selected 9). Follow-up: #691.
     - #689
     - branch ``feature/cone-store-b2``, uncommitted at the time of writing
   * - 2026-11-02
     - **Ruling 33's sub-claim refuted** ("a name collision can only
       over-select": 24 of 42 attack rows under-selected under the B2
       design) and **ruling 43** (R1 to R6, ``del`` targets as bindings);
       **ruling 44** (the skeleton retired, the store versioned).
     - #689
     - plan ``238e5103``, ``e9df7ed9``, ``549652ed``
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
     - ``7fdc7474``
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
