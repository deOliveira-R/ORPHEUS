.. _verification-cone-store:

=============================================================================
The cone store: the tests a working-tree change requires
=============================================================================

.. Machine header. Owner of: the cone of a test (the pins its verdict
   depends on) and the store that keeps two per test; why a cone is
   captured by fork with the test alone (the in-process cache leak), why a
   grandparent and one collecting parent per test module (ruling 38's
   refutation), the reading of a recording as a cone (roots, kept and
   dropped pins, the shadow presence pins, the cone-unknown reasons), the
   SQLite layout, the two-stage freshness, select's rule, the environment
   key, every read the capture silences and why, the positive controls and
   the mutation battery, the scope boundaries; since B2 (#689), the
   import-time statement as the unit of a module's import-time code (the
   two ways a run observes one, the class rules, the projection pin, the
   import-outcome pin, the quiet-section rule); and, since B3, the two
   cones of a test (the selection cone select judges, the own cone a
   second recording takes around the test's protocol), the grandparent's
   share, the store's FORMAT, and place, which homes a layer's tests in
   the layer's concept order; and, since B5, what a test's run read (the
   bytecode reader, the resolution of a load chain to the statement that
   binds its end, the channels, the B5 statement pins, the own_reads
   table), the tiered home, the apparatus nodes and the declared edges'
   c_items. Not owned here: the concept order itself
   (tests/gates/numerics/test_concept_order.py); and the pin kinds, the recorder, validation, the kept reads,
   the fork guard and the refusals (reference_cache.rst, whose sections
   this page cites); the
   dependence DAG of the suite and its rests_on edges
   (principles.rst, vv-test-suite-layout); the V&V markers (harness.rst).
   Plan of record: .claude/plans/verification_tree_refinement.md, from
   "The cone store -- design" to its end (rulings 32 to 49, ruling 38 as
   corrected, ruling 33 as amended by 43). Specification: scratch/verification_tree/cone_store/spec_b1.md;
   reviews qa_b1.md, elegance_b1.md, battery_b1.md, m1_diagnosis.md in the
   same folder; for B2, b2/gate_spec.md, b2/qa_b2_design.md (with the 42
   attack rows in b2/attacks/), b2/qa_b2_build.md and b2/elegance_b2.md;
   for B3, b3/spec_b3.md and b3/elegance_b3.md; for B5, b5/spec_b5.md,
   b5/qa/qa_b5_review.md, b5/elegance/review.md and the census of list 1,
   b4/census/l1.md.
   Written for increment B1 (landed 7fdc7474), updated for B2 (landed
   6e52b3cf), for B3 (landed 9195dcf5), and for B5 on branch
   feature/b5-placement-reads, uncommitted, 2026-11-07, over 4fad267e.

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
   python -O -m tools.cone_store place --out F   # home each stored test of a layer in its concept order

Every command takes ``--repo`` (default the working directory), ``--store``
(default ``<repo>/.cache/cone_store/cones.sqlite``), ``--package`` (the
first-party package, default ``orpheus``) and ``--memo-root`` (the traced
memo's store, default the memo's own); ``capture`` takes ``--route
{auto,fresh}``; ``place`` takes ``--layer`` (the test-id prefix of the
layer, default ``tests/gates/numerics/``), ``--code-root`` (the layer's
code, default ``orpheus/numerics``), ``--order`` (the file holding the
concept order, default ``tests/gates/numerics/test_concept_order.py``),
``--threshold`` (default 0.95) and ``--out`` (the JSON report's path;
without it the report is printed). The ids are pytest node ids or paths, as
pytest takes them. The same four verbs are functions of the package
(``tools.cone_store.capture``, ``select``, ``status``, ``place``, and
``open_store`` for reading a store). The
command must run under the interpreter and ``-O`` level the suite runs
under, since both are in the environment key.

The full suite stays the merge gate. ``select`` is a local check. The
validation run of ruling 41 (sampled mutants, each run both in full and
through the store, every selection agreeing) has run on the numerics layer
and found no under-selection (:ref:`verification-cone-store-validation`);
the full suite stays the merge gate until every layer is captured and
validated (the plan's ruling 46, the user's, 2026-11-07). The store is also the input of the census
that places tests in the verification tree (ruling 22 of the plan):
``place`` homes each stored test of one layer in that layer's concept order
(:ref:`verification-cone-store-place`).


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
* **A test has two cones** (B3). Its SELECTION cone is everything its
  process ran, the conftests' imports and its module's collection included;
  ``select`` judges it. Its OWN cone is what pytest's protocol for that one
  test ran, recorded by a second ``Recording`` started around the protocol;
  ``place`` reads it. The floor every test inherits says nothing about what
  a test exercises: placed from the selection cone, 0 of 180 sampled tests
  of ``tests/gates/numerics`` read as pure (`[M]` the B3 specification,
  2026-11-03); placed from the own cone, the Krylov test
  ``test_krylov_acceleration_with_fission_term`` is homed at ``iteration``,
  pure (:ref:`verification-cone-store-own-cone`).
* **Placement homes a test at the highest-ranked layer node of its own
  cone** (``place``) and calls it pure when every other node of its own cone is an
  ancestor of the home in the layer's import graph, read under three
  cumulative closures (top-level imports; and imports inside a def; and
  ``TYPE_CHECKING`` imports), so that each omitted kind of edge counts the
  straddles it explains (:ref:`verification-cone-store-place`). On B4's store
  of the numerics layer it placed 5 257 of 5 808 tests, 0.811 of them pure
  against a shuffled null of 0.060 to 0.064 (`[M]` 2026-11-07,
  :ref:`verification-cone-store-b4-placement`).
* **The home is tiered** (B5, ruling 49). The defs decide it when the test ran
  any layer def; else the layer bindings its run read, each load chain of
  its executed code resolved from the bytecode to the statement that binds
  its end; else the layer files it opened as data; else the test-tree
  apparatus it used. A value read never moves a def home, and the home every
  channel would give together is reported beside it as the refuted reading
  (:ref:`verification-cone-store-reads`).
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
* **Evidence.** 127 gate ids in 10 files under ``tests/gates/tools/``
  (4 of them ``slow``; `[M]` 2026-11-03, the verification matrix the docs
  build regenerates, and the ``slow`` decorators), and 91
  memo-level rows in ``tests/gates/numerics/test_traced_memo_statement_pins.py``;
  B1's 35-arm mutation battery with a positive control, every arm reddening
  its target, and a 14-arm round after review, every arm reddening its
  target and no declared-green row; B3's two batteries, each arm reddening
  its target row (:ref:`verification-cone-store-place-gates`); the memo's
  own suites give the same outcomes captured and run plainly (M1, `[M]`
  1342 s at B1).


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
       (``__pycache__`` left out of every listing, and a ``__pycache__``
       directory's own listing not recorded)
     - a name is added, removed or renamed
   * - ``ChildPin``
     - a traced-memo entry the test read from the memo's store, by its key
       and payload digest (an entry read from another store is pinned by
       its entry's own pins instead)
     - the entry no longer validates (its own manifest, checked recursively)
   * - ``DistributionPin``
     - an installed distribution whose code ran or whose file was read, by
       version and the digest of its ``RECORD``
     - another install of it is there

A row of the store is (test, environment key, pin), with the identity the
pinned file had at capture (ruling 39). A test has two such sets: its
selection cone, which this page means by "the cone" unless it says
otherwise, and its own cone, the code pytest's protocol for the test ran
(:ref:`verification-cone-store-own-cone`). The grain is the memo's: def bodies,
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

A ``ChildPin`` names no store: it is judged against the store that
``select`` validates with, the memo's default one. A test that reads an
entry from a store of its own (a ``cache_root`` under ``tmp_path``, as the
traced memo's own gates do) read an entry that store does not hold, so its
``ChildPin`` read ``Absent`` on every ``select`` and the test was always
selected (`[M]` 2026-11-03, the B3 acceptance: 70 of the 229 tests a
no-edit ``select`` listed on the numerics layer). The pin cannot simply be
dropped: the generated code ran in the generating process, which the test's
recording never sees, so without it an edit to the generator would select
nothing. A ``ChildPin`` that named its store would not help either: a
store under ``tmp_path`` is gone by the time ``select`` runs. The memo's
lookup therefore returns, with a ``Hit``, the pins it validated the entry
on, each child replaced by its own entry's pins (``Hit.pins``), and the read
is recorded with them and with its store's root (``ReadChild``).
``collect_pins`` pins a child of the collecting store by its ``ChildPin``
and a child of any other store by those pins, with no second read of the
store. A pin's source path is relative to the generating process's
``sys.path``, which was the reader's, and the entry a test put there (its
package under ``tmp_path``) is gone by collection, so each is resolved
against the ``sys.path`` recorded with the read (`[M]` 2026-11-03, qa: left
relative, the child's ``ImportOutcomePin`` re-measured its rows with no root
and read stale with no edit); a source no entry holds a file at any more is
the refusal ``SourceGone``, never a pin dropped in silence; the cone is then
unknown, since which ``sys.path`` entry held the file, and so whether it was
the repository's, cannot be told once it is gone (qa, 2026-11-03). The pins are then kept and dropped like the run's own
in the selection cone, and are kept out of the own cone unless the test's run ran them too: the generating
process ran them (``CollectedPins.through_reads``; ``cone_of`` returns the two cones as ``ReadCones``; `[M]` 2026-11-04, the numerics re-capture: with them in the own
cone, 70 of the layer's readers were homed in ``convergence`` by import-time defs of the generation, such as
``IterationBudget.__post_init__``, where their own runs home them in ``traced_memo``),
and a file only the child's generation read gets its identity recorded like
any other (qa, 2026-11-03: it was recorded as missing, which turned off
stage one for it). Gate
S23 (``tests/gates/tools/test_cone_store_selection.py``): a test reading
its own store's entry, whose generation imports a module the test never
imports, is not selected with no edit and is selected by an edit to that
module's def; the memo-side gates are in
``tests/gates/numerics/test_traced_memo_findings.py``.


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
the module grain did not depend on the names, so both gates passed. Both
pass after the repair (`[M]` 2026-11-03, the B3 gate run: the
``test_traced_memo_*.py`` files, the concept-order, axis-identity,
characteristic-block-memo, ambient-state and layer-import gates and
``tests/gates/tools/``, 1162 of 1163 passed; the one red, M1, is the
apparatus gate's premise below).

.. _verification-cone-store-apparatus:

**The recorder's apparatus is never recorded.** The quiet flag cannot
exempt the code the recorder itself runs on every event: the audit hook,
the ``os.stat`` and ``os.lstat`` probes and what they call run inside every
recording, before any ``unrecorded()`` stretch can begin, and the
``unrecorded`` and ``spawning`` stretches' own ``__enter__`` and
``__exit__`` run at its edges. So ``orpheus/numerics/_traced_memo_boot.py``
declares that code at its definition with the decorator ``_apparatus``,
into the set ``_APPARATUS`` keyed by (file, code object): ``_audit``,
``_probing`` and each probe it makes, ``_observe``, ``_path``,
``_by_importer``, ``_access``, ``_note_written``, ``_program``,
``note_child``, ``Recording.start`` and ``Recording.stop``, the
``__enter__`` and ``__exit__`` of ``unrecorded`` and ``spawning``, and
``_apparatus`` itself. A def is declared with every body nested in it (a
generator expression, a closure), because each is a code object of its own,
started on its own: ``_program``'s generator expression, the monitoring
callback ``started`` inside ``Recording.start`` and the probe inside
``_probing``. The ``PY_START`` callback returns ``DISABLE`` for
a code object in that set without recording it. The rest of the recorder's
file is recorded like any other code when a run calls it, so a test of the
recorder still depends on what it ran. The apparatus is the instrument,
pinned whole by the memo's own source (``MemoPin``) and, for the cone
store, by the environment key (:ref:`verification-cone-store-environment`),
never by a run's cone. `[M]` 2026-11-03, the same 10 sample tests captured
before and after: the recorder's defs in the own cones went from 6 to 0 in
10 of 10 tests, and in the selection cones from 7 to 2 (``Kept.__init__``
and ``_apparatus``, the module-level code that importing the recorder
runs). Gate ``tests/gates/numerics/test_traced_memo_apparatus.py``: a run
that probes a path and reads a file records no recorder def, with the
probe's observation as the activation leg; a run that calls the recorder's
``digest`` records it; every hook a recording installs is read off the
installed objects and found declared (emptying ``_APPARATUS`` reddens the
first and third, declaring every code object of the file apparatus the
second); a nested recording run over every hook path records no apparatus
in either recording; the same nested run with the audit hook made
traceable records none either (its first red: ``_program``'s generator
expression, while ``_apparatus`` declared a def's own code object and not
the bodies nested in it); and 16 arm cases, one per declared def, each
remove one def with its nest and require the nested run to record that
def, so every declaration is reached by some hook path. The traceable run
and the arm cases each run in a fresh interpreter, for two reasons. A
recording installs the ``os.stat`` probes, and so calls ``_probing``, only
when no recording is active, so the first recording of a process is the
only one that reaches ``_probing`` and ``_apparatus``; in-process under the
cone store's capture, whose child records every test, those two arm cases
read red while green in a plain run (`[M]` 2026-11-03, gate M1). And a
traceable audit hook set in-process would reach the recording around the
test. The audit hook and the defs only
it calls (``_access``, ``_note_written``, ``_program``) are a special case:
CPython calls an audit hook with ``sys.monitoring`` paused unless the hook
sets ``__cantrace__``, so today they are never started under a recording
at all. They are declared anyway, and the traceable run and their arm
cases set ``__cantrace__`` on ``_audit`` to witness the declarations. Editing the recorder changes the memo's
source digest, so every traced-memo entry goes stale once and regenerates
with the same values.

.. dropdown:: First got wrong: "the recorder's hooks are the capture's instrument, counted apart by placement"
   :color: muted

   The recorder recorded its own hooks: with any recording active, the
   audit hook and the stat probe ran inside each test's protocol, so the
   same six recorder defs (``_probing``, ``_path``, ``_observe``,
   ``_by_importer``, ``unrecorded.__enter__``, ``unrecorded.__exit__``) were
   in every cone the store read. The first ``place`` treated the recorder's
   file as an instrument region, its defs counted and never placed, with a
   fallback that homed a test at the recorder when its own cone held no
   other layer def, meant for a test of the recorder itself. The fixture
   could not show what that fallback did: its stand-in recorder was called
   explicitly, and the real recorder's file lay outside the fixture's layer.
   The archivist, documenting B3, read the placement of the orchestrator's
   10-test sample store: `[M]` 2026-11-03, 10 of 10 own cones held the six
   defs, and 6 of the 10 tests were homed at ``_traced_memo_boot``: the 5
   tests of ``test_concept_order.py``, which read the layer's sources and
   run none of its code, and
   ``test_krylov_acceleration_requires_a_stated_preconditioner``, whose
   refusal is the interpreter's argument check. The root cause was the
   recorder violating this section's rule, so the repair is at the
   recorder (the apparatus above); the instrument region, the fallback, its
   count and its gate PL12 were retired, and ``_traced_memo_boot`` is an
   ordinary ranked node. After the repair the 6 tests are unplaced with "no
   layer def" and the 4 placed tests read as before.

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
   as ``Corrupt`` and regenerates; the cone store's rows carried
   ``FORMAT = 2`` (3 since B3, 4 since B5) and a store of another form is emptied when
   opened (:ref:`verification-cone-store-store`). An
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
   child that writes the grandparent's cone, the floor it attributes to
   every test (stored once, as a share:
   :ref:`verification-cone-store-own-cone`), then forks one collecting
   parent per test module, waiting for each.
2. **A collecting parent** (``Collector``) is that fork continuing as an
   ordinary pytest session over its module's arguments. It collects the
   module and, for each test, forks a child.
3. **A child** runs pytest's own protocol for the one test (setup, call,
   teardown) under a second recording of its own, stops both recordings,
   reads each as a cone (:ref:`verification-cone-store-reading`) and writes
   the two, with its route and outcome, to the capture's results directory
   (:ref:`verification-cone-store-own-cone`).

The recording is inherited across both forks and only grows, so a child's
selection cone is the grandparent's (the conftests' imports) plus its module's
collection plus its own run. Nothing re-arms the monitor: a code object the
recorder disabled after its first start is already in the set every later
child inherits. The own cone is a second recording, not a re-arming of this
one (:ref:`verification-cone-store-own-cone`).

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
still recorded. The recorder's own apparatus is never recorded at all
(:ref:`the recorder's apparatus <verification-cone-store-apparatus>`).

**The capture's own code is in no cone.** The worker module runs in every
capture and in no full run, so it is no test's dependency. `[M]` 2026-11-02
(a capture of ``tests/gates/test_layer_imports.py``, 533 tests): six of its
defs (the hooks pytest calls inside a recorded stretch, and the wrappers
that switch the recording) were pinned in every cone, and so were the
import-time statements of ``cone.py``, ``pins.py``, ``store.py`` and
``worker.py``, so an edit to the store's own code selected the whole suite
(qa's and the elegance-enforcer's final pass on B2, F5). The worker removes
its own module's code objects from each of its recordings before it reads
it (the selection one and, since B3, the own one), and imports the modules
that read a cone only after both recordings stopped. A
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

.. _verification-cone-store-own-cone:

Two cones per test: the selection cone and the own cone
-------------------------------------------------------

The cone ``select`` judges is the SELECTION cone: everything the test's
process ran, which by rulings 38 and 42 includes what the grandparent ran
while it imported the conftests and what the module's collecting parent ran
while it collected. That attribution is right for selection, since a full
run runs that code before every test. It is wrong for any question about
what the test itself exercises, because the floor is the same for every
test. On a sample of 180 tests in five files of ``tests/gates/numerics``,
139 ``DefPin`` rows are common to all 180 selection cones, of 142 to 159
per test, and a placement read from the selection cone finds 0 of 180 tests
pure, 139 of them homed at ``convergence``, the module of
``IterationBudget.__post_init__``, which the floor runs (`[M]` the B3
specification, ``scratch/verification_tree/cone_store/b3/spec_b3.md``,
its measurements M8 and M9, 2026-11-03).

So each child keeps two cones, stored as two relations over the same pins:

* the **selection cone**, from the recording the child inherited (the
  grandparent's, continued through both forks). It is the cone B2 stored,
  row for row, and ``select`` reads it and nothing else, so selection is
  unchanged;
* the **own cone**, from a second ``Recording`` that the child starts
  immediately before pytest's ``pytest_runtest_protocol`` for its test and
  stops when the protocol returns (setup, call and teardown;
  ``Recorded.run`` in ``tools/cone_store/worker.py``). ``place`` reads it
  and nothing else.

**Why a second recording sees code the first has already seen.** The memo's
recorder is a ``sys.monitoring`` ``PY_START`` callback: it records a code
object the first time the object starts, then returns
``sys.monitoring.DISABLE``, which turns that event off for that code object
and for that monitoring tool only. Each ``Recording`` takes a tool id of its
own when it starts (the first free id among 2, 3, 4 and 5;
``Recording.start`` in ``orpheus/numerics/_traced_memo_boot.py``), so on a
code object the outer recording has disabled, the inner recording's event
is still enabled. A def the conftests ran at import, and that the test runs
again, is therefore reported to the own recording when the test starts it;
a def only the conftests ran is not in the own cone. Gates PH1 and PH2 pin
the two cases on a fixture whose root conftest runs ``conftest_helper`` at
import: ``test_p`` calls it again, and its own cone holds it; ``test_b``
never does, and its own cone lacks it while its selection cone holds it.
With one recorder serving both cones, the own cone would lack every def the
floor had started first (PH1's first red). Stopping a recording frees its
tool id and calls ``sys.monitoring.restart_events()``, which re-enables
every disabled event; the outer recording then reports again code it
already holds, which adds nothing, since its record is keyed by file and
code object. A test that nests a ``trace_call`` takes a third id; a fifth
recording active at once would find no free id and raise.

**The reading.** Both recordings are read after both have stopped (the
selection one first), and the modules that read a cone are imported only
then, so their import-time code is in neither cone. The worker removes its
own module's code objects from each recording before reading it (above,
"The capture's own code is in no cone"), and each recording is read by the
same ``cone_of`` and kept by the same rules
(:ref:`verification-cone-store-reading`), which returns both cones of the
recording (``ReadCones``): the own cone is the ``ran`` cone of the own
recording, without the pins the run holds only through a memo entry read
from another store, which another process ran; a pin the run ran itself
stays, whatever entry also holds it. Two cones of one child name a def
by the same ``DefPin`` row, but not a module by the same ``StatementsPin``
or ``ImportOutcomePin`` row: a statement pin is a projection on the entry
names of the run that recorded it, the own run's entry names are fewer, and
its import set is shorter. So containment between the cones holds by NODE
(the pin's kind, path, qualname and ordinal), never by row: gate PH3
asserts, for every fixture test, that the nodes of its own cone and of the
grandparent's share lie in the nodes of its selection cone, and that an
edit to ``conftest_helper``'s body selects ``test_b``, which reaches the def
only through the conftest (its first red: ``select`` judging the own
cones).

**Why two relations, and not a phase per row.** The B3 specification ruled
a phase on each row (``inherited`` or ``own``), with ``select`` reading the
union. Two relations were built instead. The selection relation stays the
relation B2 wrote: the B2 tools of ``83db308f`` and the B3 tools captured
the same 7 tests (``test_concept_order.py`` and two of
``test_iteration.py``) with equal selection rows for 7 of 7, 1 860 to 1 943
rows each (the one-off comparison recorded in PH3's docstring, 2026-11-03).
And the own cone is a sub-recording of the selection cone, not a partition
of it: a def both the floor and the test ran belongs to both, which one
phase per row cannot say (the elegance review of B3, its rejected attack
"a phase column").

**Why not subtract the floor.** Removing from the selection cone the pins
every test shares, or the grandparent's stored share, loses every def the
test shares with the conftests: on the sample at least
``IterationBudget.__post_init__``, which a Krylov test constructs itself.
The B3 specification's prototype placed tests that way (its ``FLOOR=1``
reading); it is refuted for the question "what did the test run", and gate
PL6 is its red (a test whose every def is also in the floor is left with no
layer def).

**The grandparent's share.** The grandparent's cone, the floor, is written
by its short child and stored once, under the scope ``"<grandparent>"``
(``store.GRANDPARENT``) of the ``shares`` table, never repeated per test;
``status`` derives the grandparent's pin count from it. A capture that
produced no grandparent cone (the fresh route, or a grandparent that died)
clears the stored share rather than keep another capture's. Gate PH4 asserts
that the share holds ``conftest_helper`` and lies, by node, in every test's
selection cone (its first red: the share not stored).

**The cost.** `[M]` 2026-11-03, a profile of one capture child:
``collect_pins`` takes 1.8 to 1.9 s per child over 4 827 code objects, most
of them the inherited floor, and reading the own recording adds 0.10 s.
`[R]` At that rate the 5 780 tests of ``tests/gates/numerics`` take about
2.6 h to capture. Every child re-reads the floor its collecting parent has
already read; amortising that read per collecting parent is #696.


.. _verification-cone-store-reading:

Reading a recording as a cone
=============================

``cone_of`` (``tools/cone_store/cone.py``) turns a stopped recording into a
cone in three steps. A capture child reads both of its recordings this way
(:ref:`verification-cone-store-own-cone`).

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
       when it is a content pin (``DefPin``, ``StatementsPin``, ``DataPin``);
       the three relations below share it
   * - ``files``
     - one row per judged file: its identity at capture and whether that
       identity is trusted
   * - ``cones``
     - (test id, pin id): each test's selection cone, the one ``select``
       judges
   * - ``own_cones``
     - (test id, pin id): each test's own cone, what pytest's protocol for
       the test ran; read only by ``place``
   * - ``own_reads``
     - (test id, path, name, channel): what each test's own run read (B5),
       one row per binding by the repository path of the module whose
       statement binds it, the bound name (empty for a read of the module
       object) and its channel; read only by ``place``
       (:ref:`verification-cone-store-reads`)
   * - ``shares``
     - (scope, pin id): a cone common to many tests, stored once; the
       grandparent's under the scope ``"<grandparent>"``
   * - ``meta``
     - the store's ``format`` and the number of tests the last capture
       pruned
   * - ``rows`` (a view)
     - the flat reading of the selection cones, one row per test and pin,
       with the file's identity

**A file's identity is a fact about the file**, stored once and shared by
every pin on it. A capture that reads a file whose stored identity differs,
or was not settled, stores it untrusted, so the next ``select`` judges every
pin on that file (stage two); where both cones of a test name a file, the
state its selection cone recorded is kept (one capture, one identity per
file). **A capture replaces a test's rows whole**, in both relations, and
replaces the grandparent's share. It
forgets tests that no longer exist: those of a re-captured module that the
module no longer holds, and, after a capture of everything, those of a module
that is gone. A capture of named node ids prunes nothing else, and that is
safe: ``select`` iterates collected ids only, and a re-added test's old row
is judged by its pins' digests against the current tree. After a capture, a
pin no relation references (``cones``, ``own_cones`` or ``shares``) is
deleted, and so is a file no remaining pin judges.

**The store's form is versioned.** ``tools/cone_store/store.py`` declares
``FORMAT = 4`` (since B5; 3 was B3's and 2 B2's) and keeps it in ``meta``. A store
opened by ``ConeStore`` with another ``format``, or with none (every store
B1 wrote), has its tables dropped and recreated in one transaction, and
every test is then new: selected, and captured again. A B1 store holds
``ModulePin`` rows, a kind the memo no longer decodes (`[M]` the elegance
review of B2: ``Manifest.from_json`` of a ``["ModulePin", …]`` row raises
``KeyError: 'ModulePin'``, which ``select`` did not catch); a B2 store holds
no own cones, so a placement read from it would find every own cone empty;
a B3 store holds no ``own_reads`` rows, so a placement read from it would
find no test's reads.
Emptying is safe in the direction that matters for ``select``: it selects
everything. Gate SG5 writes a store in B1's form and requires ``select`` to
list every collected test (dropping the format check reddens it with that
``KeyError``); PH5 writes a store of format 2 and requires it emptied when
opened, with ``own_cones`` and ``shares`` created; P8 marks a captured store
with form 3 and requires its reads emptied when it is opened, and
``FORMAT`` to be 4.

``place`` only reads a store, so it must not empty one. It reads the format
first with ``stored_format(path)``, which opens the database read-only
(SQLite's ``mode=ro``) and returns the stored format, or ``None`` when
there is no store or no ``meta`` table, and it refuses a store of any
other form by name, saying what that form lacks (``_LACKS`` in
``placement.py``): a store of form 2 "holds no own cones (written before
B3)", one of form 3 "holds no own reads (written before B5)", and a store
of any other form is said to lack "the rows this placement reads"; each
message ends "so there is nothing to place; capture again". Gate PL11 copies a captured
store, sets its format to 2, and requires the refusal naming the own cones.


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
environment key, distinct pins and selection rows, the grandparent's pin
count (the size of its stored share), the
last capture's pruned count, and each unknown test with its reasons.


.. _verification-cone-store-place:

Placing a layer's tests: ``place``
==================================

``python -O -m tools.cone_store place`` reads the own cones of a store, and
since B5 the own reads (:ref:`verification-cone-store-reads`), and
homes each stored test of one layer in that layer's concept order: the
census step of ruling 22 of the plan. The code is
``tools/cone_store/placement.py``. The command loads the concept order from
a file by path (``--order``, default
``tests/gates/numerics/test_concept_order.py``), reading two names from
it: ``TIERS``, the order, and ``import_kinds``, the import graph. The report
is one JSON document, written to ``--out`` or printed; with ``--out`` the
summary is also printed, with the straddler list as its count.

The report places; it decides nothing. Keeping, strengthening or retiring a
test is phases 2 and 3 of the plan and ruling 21's mutant method, and the
duplicate candidates below carry no verdict. Numerics has no layer below it
whose claims are placed, so its tests are placed against the module order
only, and its proposed supports run between numerics tests.

Nodes, regions and the order
----------------------------

A layer's code is cut into **nodes** at the concept order's grain: a module
directly under the layer's root (``orpheus/numerics``) is a node, and a
package directly under it (``basis``, ``quadrature``, ``spaces``) is one
node; the layer's own ``__init__.py`` is no node. The order is a sequence of
tiers, upstream first; a node's rank is its tier's index, and the nodes of
one tier share a rank. The recorder, ``_traced_memo_boot``, is an ordinary
node, and a test's own cone holds a recorder def only when the test ran
one (:ref:`the recorder's apparatus <verification-cone-store-apparatus>`). `[M]` 2026-11-03 (a probe reading ``TIERS`` and
``import_kinds()``): the numerics order has 50 nodes in 25 tiers.

Every ``DefPin`` of a test's own cone falls in one **region**, a sum type
(``Region`` in ``placement.py``) consumed by one exhaustive ``match``
(``own_defs``):

.. list-table:: Where an own def is
   :header-rows: 1
   :widths: 18 42 40

   * - Region
     - Which defs
     - What ``place`` does with them
   * - ``InLayer(node)``
     - a def in a file of the layer, under its node
     - the test's layer nodes, which it is placed by
   * - ``Upstream(root)``
     - a def under a layer below this one (``orpheus/geometry`` and
       ``orpheus/values`` for numerics)
     - counted per module under ``upstream_layers`` (``orpheus/geometry/<node>``;
       the layer's root itself for its own ``__init__.py``); never a
       straddle
   * - ``Above(package)``
     - a def under any other directory of the first-party package
       (``orpheus/data``, ``orpheus/transport``, a method)
     - counted per directory under ``above_layer``: a test of this layer
       reaching a layer above it, a finding of its own
   * - ``Helper()``
     - a def under ``tests/``: a helper, a fixture, a test-tree def another
       test file exports
     - kept per file as the apparatus the test ran, its items (B5), and
       counted under ``helpers``; the test's own def, the code nested in it
       and the pytest hooks are not apparatus and are not counted (in B3
       every def under ``tests/`` was counted)
   * - ``Elsewhere()``
     - first-party code in no region above: the layer's own ``__init__.py``,
       a file outside the package (``tools/``)
     - counted under ``elsewhere``

The regions sort defs. B3's placement read defs only, and a data, presence
or listing pin was a read and never a node, so the tests of
``test_concept_order.py``, which read the layer's sources and run none of
its code, were unplaced. Since B5 the def tier is the first of four: a test
that ran no layer def is homed by the layer bindings its run read, else by
the layer files it opened as data, else by its apparatus
(:ref:`verification-cone-store-reads`). The rest of this section is the def
tier, which is B3's unchanged.

Home, tie and purity
--------------------

For a test ``t`` placed by its defs, ``own(t)`` is the set of layer nodes
that hold a def of its own cone. A test placed by its reads takes its home
and tie from its statement or data nodes by the same rule.

* The **home** is the node of ``own(t)`` of highest rank. When two or more
  nodes share that rank, the test has a **tie**: the report gives the tied
  nodes and no home, and the tie is never broken (by name, by def count, by
  order of appearance). A tied test is pure under no closure, and its tied
  nodes are listed among its foreign nodes with the kind ``tie``.
* The test is **pure** when every other node of ``own(t)`` is an ancestor
  of the home in the layer's import graph: the home's module imports it,
  directly or through other nodes.

The closure is the import graph and never the rank order. Under the rank
order every node below the home would count as upstream, purity would mean
only "no tie", and the shuffled null (below) would read the same as the
placement.

A **straddler** is a placed test that is not pure under the runtime closure.
Each of its **foreign** nodes (in ``own(t)``, not the home, not a runtime
ancestor of it) carries the kind of the first closure that explains it:
``lazy`` (an ancestor once imports inside a def count), ``typing`` (an
ancestor once ``TYPE_CHECKING`` imports count too), or ``none`` (an
ancestor under no closure: the test composes the two nodes itself).

Three closures of the import graph
----------------------------------

Each edge of the import graph carries the kinds of the import statements
that make it, and ``place`` reads purity under three cumulative closures,
side by side:

.. list-table::
   :header-rows: 1
   :widths: 14 14 72

   * - Closure
     - Report key
     - Edges admitted
   * - ``runtime``
     - ``pure_rt``
     - top-level imports, wherever they sit in the module outside a def:
       the imports the interpreter runs when the importer is imported
   * - ``lazy``
     - ``pure_lazy``
     - and imports inside a def or an ``async def``, run when the def is
       called
   * - ``typing``
     - ``pure_typ``
     - and imports under ``if TYPE_CHECKING:``, which never run; every import
       in such a block is ``typing``, even inside a def

The kinds are kept apart because they are different dependences. A runtime
import holds whenever the home's module is loaded. A lazy import holds only
when the importing def runs. A ``TYPE_CHECKING`` import is a name the type
checker reads and the interpreter never runs, so a test that reaches a node
only along such an edge ran that node's code because the test built the
object itself, not because the home ran it. The three columns side by side
count the straddles each omitted kind of edge explains: a test pure under
``typing`` and not under ``runtime`` straddles only along edges the
interpreter does not follow at import. `[M]` 2026-11-03 (the same probe):
the numerics graph has 118 edges between distinct nodes; 90 carry a
top-level import, the lazy closure admits 12 more, and the typing closure
the last 16.

**One walk of the import graph.** The graph is ``import_kinds`` in
``tests/gates/numerics/test_concept_order.py``, the concept-order gate's
own walk: for each pair (imported node, importing node) of distinct nodes,
the set of kinds that make the edge. The gate reads its edges
(``import_edges`` is the set of its keys); ``place`` receives the kinds as
an argument. ``tools/cone_store`` parses no source: gate R7 forbids
``ast.parse`` there, so no second walker can be written in it; the walkers
the tree still holds are #695 (one import-graph instrument). The walk
visits each statement with the kind of its scope: ``runtime`` at module
level, ``lazy`` inside a ``def`` or an ``async def`` (a lambda cannot hold
an import statement), and
``typing`` in the body of an ``if`` whose test is exactly
``TYPE_CHECKING`` or ``<name>.TYPE_CHECKING``; the ``else`` of such an
``if`` keeps the enclosing kind. The predicate is exact, not a substring
test, so ``if not TYPE_CHECKING:``, which guards the runtime branch, is
read as an ordinary ``if``; a substring test reads its body as ``typing``
(the elegance review of B3, C3; no such spelling is in ``orpheus/``
today). The gate ``test_each_import_edge_carries_its_kind`` plants one
import of each kind, a ``TYPE_CHECKING`` import inside a def and an edge
made by two kinds, and requires the walk to return exactly those kinds;
``test_import_kinds_reads_the_branch_the_interpreter_runs`` plants an
import in the ``else`` of ``if TYPE_CHECKING:`` (``runtime``), one under
``if not TYPE_CHECKING:`` (``runtime``) and one in an ``async def``
(``lazy``), each reddened by its textual mutant of the walk.

**An unknown kind is refused.** The kinds are spelled twice: by the walk
(its ``KINDS``) and by ``placement.CLOSURES``, which maps each kind to the
closure that first admits it. A kind the walk gained without a closure here
would read as ``none``, silently, so ``place`` refuses it before reading
the store ("import kinds this placement has no closure for: [...]"). Gate
PL14 hands it the kind ``conditional`` and requires the refusal.

**Ancestors are reachability, so a cycle closes exactly.**
``ancestors(kinds, nodes, admitted)`` builds the graph of the admitted
edges, importer to imported, as a sparse matrix, and runs one breadth-first
order per node (``scipy.sparse.csgraph.breadth_first_order``): a node's
ancestors are the nodes reachable from it, itself excluded, whatever order
the nodes are visited in. The numerics graph has a cycle:
``traced_memo`` imports ``_traced_memo_boot`` at top level, and the boot
imports ``traced_memo`` inside a def so that the memo can boot itself (the
concept-order gate's one ``EXEMPT`` edge). `[M]` 2026-11-03 (the probe):
under the runtime closure the boot is an ancestor of ``traced_memo`` and
not the reverse; under the lazy and the typing closures each is an ancestor
of the other.

.. dropdown:: First got wrong: "a memoised depth-first walk gives each node's ancestors"
   :color: muted

   The first ``ancestors`` was a depth-first walk that memoised each node's
   closure when it finished. On a cycle it stores a closure computed while
   part of the cycle was still being visited, so its answer depended on the
   order the nodes were visited in. On the live graph, ``traced_memo``'s
   typing closure lacked ``_traced_memo_boot`` while its runtime closure
   held it, so a support onto a test homed at the boot would have been
   listed as against the order (the elegance review of B3, V1). Reachability
   over a sparse graph replaced it. Gate PL13 runs ``ancestors`` on three
   nodes with a 2-cycle in three visiting orders, and on the live pair;
   `[M]` the memoised walk reddens it in 3 of 3 orders.

Unplaced tests
--------------

A test is **unplaced**, and kept in the report with its reason, when (checked
in this order):

* its cone is unknown ("cone unknown: " and the stored reasons): its record
  stops where the recording lost sight of the run;
* it did not pass at capture ("not passed: " and the outcome, ``failed``,
  ``skipped``, ``crashed`` or ``uncaptured``): its cone stops at the
  failure;
* it read nothing of the layer or of the apparatus ("no layer read"): no
  layer def ran, no layer binding was read on a node channel, no layer file
  was opened as data, and no test-tree item was used. B3's reason was "no
  layer def", met by every test that ran no layer def; since B5 most such
  tests are placed by their reads.

The summary counts them by the text before the colon. An unplaced test keeps
its region counts and carries no home.

The shuffled null
-----------------

A placement that a random labelling matches carries no information (S13 of
the plan). The null takes every distinct (node, def) pair of the placed
tests' own layer defs, shuffles the node labels over those pairs with
``random.Random(seed)`` for each of three fixed seeds (1, 2 and 3), and
reads the runtime purity of every placed test again under the shuffled
labels. The summary gives the measured pure fraction beside the null's
fraction per seed. The seeds are fixed so that the report is byte-identical
across runs (gate PL10 compares two runs' JSON with sorted keys). A
placement that separates reads well above its null: on the specification's
180-test sample, placed with the floor subtracted (the prototype, not this
code), 0.79 against 0.05 to 0.06 (its measurement M10).

Proposed supports
-----------------

``place`` also proposes ``executes`` supports between tests (S17 of the
plan): a support ``U → L`` says that the upper test ``U`` executes nearly
everything the lower test ``L`` executes, so that ``L`` is a candidate rung
under ``U``.

* **The unit is the test function.** The cases of a parametrised test are
  one claim (ruling 4 of the plan), so a test id is cut at its first ``[``
  and the own defs of its cases are united.
* **The defs are the package's.** ``D(T)`` is the set of a function's own
  defs in the layer, in a layer below it or in a layer above it, each as
  (path, qualname). Test helpers are left out: a helper is unique to its
  file, so containment across files would be impossible by construction;
  so is code ``Elsewhere``.
* **Containment is fractional:** ``c(U ⊇ L) = |D(L) ∩ D(U)| / |D(L)|``. A
  support is proposed when ``c ≥ threshold`` (inclusive; 0.95 by default),
  ``|D(L)| ≥ 3`` (``min_defs``: a cone of one or two defs is contained by
  nearly everything), and the two cones are not equal.
* **The lower's subject is contained at the same threshold** (ruling 50 of
  the plan). The subject ``S(L)`` is the lower function's own defs at its
  home nodes (each tied node for a tie; each case's home, its cases united,
  for a function whose cases home apart), and the support also needs
  ``c_home(U ⊇ L) = |S(L) ∩ D(U)| / |S(L)| ≥ threshold``. The reason: ``D``
  holds the layers below too, so a lower whose cone is one numerics def and
  twenty ``geometry`` defs is contained at 0.95 by an upper that ran the
  twenty and not the one. The 5 % the whole cone may miss would then fall
  entirely on the code the lower is about, and a defect there cannot reach
  the upper, so the upper does not rest on the lower. ``c`` keeps its meaning
  (the whole cone); a support record carries both.
* **Equal cones are duplicate candidates** when they hold at least
  ``min_defs`` defs, listed with their def count and no verdict, and never a
  support: ruling 21 decides on kill sets, not on cones. Two smaller equal
  cones are below the size a containment is read at, and are not listed.
* **Mutual containment is a near-duplicate group.** Fractional containment
  can cycle (``U`` holds 95 % of ``L`` and ``L`` 95 % of ``U``). The
  functions that contain one another at the threshold form strongly
  connected groups (``scipy.sparse.csgraph.connected_components``, strong
  connection); a group of two or more is listed as near-duplicates, with no
  verdict, and no support is proposed inside it.
* **The proposals are the transitive reduction over the condensation.**
  Each group becomes one node, and a support between two groups is kept
  only when no other successor of the upper group reaches the lower group.
* **Supports against the order are listed apart.** A kept support whose
  lower function's home is neither the upper function's home nor an
  ancestor of it under the widest (``typing``) closure goes to
  ``supports_against_order``, for the orchestrator's reading (S18 of the
  plan), not to ``supports``.

The report also counts the containing pairs, duplicates excluded, at the
thresholds 1.0, 0.95, 0.9 and 0.8, so that the threshold can be judged on
the population it is applied to. Why 0.95: on the specification's sample
of 77 declared ``rests_on`` edges between numerics tests, 26 held strict
containment and 39 held 0.95 or more, followed by a tail down to 0.002 that
is reads and supports edges, not executes (S17); 0.95 is that
distribution's knee. The sample's edges are not the layer's declared
population, and the tail is not read as a second kind of edge: a declared
``rests_on`` is an execution support (ruling 51 of the campaign plan), so
an edge in the tail is a defect of its marker, and over the layer's 443
declarations 170 hold at 0.95 ("What a declared edge asserts", below). On the
specification's 180-test sample (with the floor subtracted), lowering it to
0.8 triples the containing pairs, 1 227 to 3 706, about 20 per test (its
measurement M12).

**The declared edges are measured on the same cones.** So that
the threshold is judged on the layer's whole declared population, ``place``
takes the declared edges (``declared``: each an upper test id and a target
id, as pytest resolved the marker) and reports ``c`` for each on the same
function cones the supports use. A case, as upper or target, collapses to
its function; a class or file target expands to every function under it;
an edge from a function to itself (a class or file target holding the
upper) is skipped. A function is measured when its test is judged: it
passed with a known cone, whatever homes it (a failed or unknown cone stops
early, so its containment means nothing). The report adds ``declared_rests_on``
(each ``{upper, lower, target, c, c_home, c_items}``; ``c_items`` and the
rows with ``c`` = ``None`` are B5's, :ref:`verification-cone-store-reads`;
``c_home`` is the subject's containment of ruling 50, ``None`` for a lower
not placed by its defs, and filters nothing: a declared edge is measured,
not proposed),
``declared_counts`` (cumulative over ``c`` at 1.0, 0.95, 0.9 and 0.8, and
``below``), ``declared_items_only`` and ``declared_unmeasured``, counted by
reason: ``upper not judged``, ``target not judged`` (no judged function
under it) and ``lower ran no def and read nothing`` (B3's reason was
``lower ran no def``). The command collects the edges through
the suite's one resolver, ``tests/_harness/rests_on_manifest.py``, loaded
as a plugin into a fresh ``pytest --collect-only`` over the layer. Gate
PL15 declares five edges on the fixture (a 0.9 edge, an edge onto a case,
two edges from the two cases of one function, a target no stored test has)
and reddens under no case collapse and under strict thresholds.

**What a declared edge asserts.** ``@pytest.mark.rests_on(L)`` on a test
``U`` asserts an execution support (ruling 51 of the campaign plan): ``U``'s
run contains ``L``'s subject, so a defect that reddens ``L`` reaches ``U``,
and a red ``L`` explains a red ``U``. That is the localisation
:doc:`principles` describes ("A failure is then localised by the lowest red
rung"). A reading order of a ladder, in which ``L``'s law is a step of
``U``'s argument while ``U`` runs other code, is not a ``rests_on`` edge:
nothing measures it, and a defect in the code ``U`` skipped cannot redden
``U``. A declared edge and a proposed support are one relation with one
predicate (ruling 53): ``c`` and ``c_home`` both at the threshold. The
defs a lower runs to build its fixture or reference stay in its cone, since
a defect there reddens ``L`` and not ``U``. A lower that ran no layer def
(it reads declared structure) is a structural premise, measured by
``c_items``. A class or file target holds for an upper when the
upper contains some lower under it (ruling 52): a marker written on a class
or in a file's ``pytestmark`` is shorthand, and each test it covers must
hold on its own. A declaration that fails is a defect of its marker,
repaired by re-pointing, narrowing or removing it.

``place`` reports each declaration's verdict in ``declared_verdicts``
(``{upper, target, holds, held_by}``, one per upper function and target as
written, for every declaration with a measured lower; ``held_by`` lists the
lowers that hold it), and each row of ``declared_rests_on`` its own
``holds``. A declaration holds when some lower under its target holds. The
containment predicate is one function, ``_supported``, called by the
proposed supports and by the declared verdicts alike: a lower with a layer
subject holds when ``c`` and ``c_home`` both reach the report's threshold;
a lower with none ran no layer def (a test that runs one is placed by its
defs, so it has a subject) and holds when ``c_items`` reaches it. Two
filters apply to proposals only: the lower's cone holds at least
``min_defs`` defs, and the two cones are unequal (an equal pair is a
duplicate). A declaration whose target holds only the upper is counted in
``declared_unmeasured`` as ``target holds only the upper``. Gates PL26 to PL31
pin it on the fixture's ``_VERDICT`` layer: a lower whose subject is held
and whose cone is not (PL26), a class target held by one of its two lowers,
``held_by`` exact (PL27), a structural lower held and missed by its read
items (PL28), the report's threshold rather than a fixed one (PL29), a
declared pair the proposals refuse on ``c_home`` (PL30), and a target
holding the upper alone (PL31).

`[M]` 2026-11-08 (census list 3,
``scratch/verification_tree/cone_store/b5/census/l3.md``, on the numerics
store at ``af5b78c8`` under ruling 50): the layer holds 443 declarations
with a measured lower (an upper function and a target as written), 976
rows; 54 more have an upper (29) or a target (25) that no judged function
answers. Where the marker is
written, read from the test files' syntax trees: on the function 199, on a
class 26, in the file's ``pytestmark`` 218. Under rulings 51 to 53, 170
declarations hold and 273 fail; 6 of the 443 are measured by read items
only (2 hold, 4 fail), and a structural lower holds by ``c_items`` also
inside a class target that mixes it with lowers that ran defs. Of the
273, 43 are ``traced_memo`` declarations that miss 1 to 4 lookup-path defs
of subjects of 17 to 37 defs, or of the cone: the template tests generate in a child forked
from the template process, which the own recording does not see, so these
wait on the capture of the child
(`#701 <https://github.com/deOliveira-R/ORPHEUS/issues/701>`__). All 62 rows
against the concept order come from markers written on a file (61) or a
class (1), inherited by tests homed at a lower node than the target.

.. dropdown:: First got wrong: "a one-hop transitive reduction is enough"
   :color: muted

   The first reduction dropped a support ``U → L`` whenever some ``M`` had
   ``U → M`` and ``M → L``, over the relation as it stood. Fractional
   containment can cycle, and on a cycle each member is a middle for the
   others: with ``A`` and ``B`` containing each other at the threshold, an
   upper ``S`` containing both lost its edge to ``A`` (through ``B``) and its
   edge to ``B`` (through ``A``), and was left with no support at all. The
   reduction now runs over the condensation, each strongly connected group
   one node, and the groups are listed as near-duplicates. Gate PL9's
   fixture keeps its two uppers out of each other's containment and names
   the cycle case in its docstring.

The report
----------

Per test (``tests``): ``home``, ``tie``, ``pure_rt``, ``pure_lazy``,
``pure_typ`` and ``foreign`` (node to kind) when placed; and always
``upstream_layers`` (module of a layer below to def count) and
``above_layer`` (directory to def count),
``helpers``, ``elsewhere`` and ``unplaced`` (``null`` or the reason). The
``summary``: the test and placed counts, the unplaced counts per reason, the
pure counts per closure, the runtime pure fraction and the null's per seed,
the ties, the foreign kinds, the homes per node (a tie as its nodes joined
by ``+``), the number of tests with an above-layer def, and the straddlers,
sorted by their foreign kinds. Beside them: ``supports`` and
``supports_against_order`` (each ``{upper, lower, c, c_home}``), ``duplicates``
(``{a, b, defs}``), ``near_duplicates``, ``support_counts`` and the
declared-edge keys above.

B5 adds, per test: ``placed_by`` (``def``, ``statement``, ``data``,
``layer``, ``package`` or ``apparatus``), ``via``, ``items``,
``home_union``, ``statement_straddle`` (for a def- or statement-placed
test), ``unresolved`` and ``apparatus_items``; to the summary:
``placed_by`` (counts), ``homes_union``, ``union_moves``,
``statement_straddlers``, ``above_layer_statement``,
``pure_fraction_statement`` and ``null_pure_fraction_statement``; and to
the result, ``apparatus_supports`` (:ref:`verification-cone-store-reads`).
The summary's ``placed`` still counts the tests placed by their defs only;
``placed_by`` counts every tier. ``place`` on a path with no store raises ("no
store at …"). The report
is plain JSON, byte-identical across runs when dumped with sorted keys.

The package function ``tools.cone_store.place`` takes the store, the
test-id prefix, the code root, the order and the kinds, and places against
the numerics layer's surroundings (``placement.NUMERICS`` with the given
root: the layers below are ``orpheus/geometry`` and ``orpheus/values``);
``placement.place`` takes a ``Layer`` for another layer or a fixture.

What the order reads on the first witnesses
-------------------------------------------

`[M]` 2026-11-03, a store of 10 tests of ``tests/gates/numerics`` captured
by the B3 tools and placed by this code, after the apparatus repair: 4
placed, 6 unplaced with "no layer def", B3's reason (the 5 tests of
``test_concept_order.py``, which read the layer's sources and run none of
its code, and ``test_krylov_acceleration_requires_a_stated_preconditioner``,
whose refusal is the interpreter's argument check). The placed four:

* ``test_iteration.py::test_krylov_acceleration_with_fission_term``: home
  ``iteration``, pure under all three closures;
* ``test_iteration.py::test_keigenvalue_recovers_dominant_eigenvalue``:
  home ``green_operator``; ``eigenvalue`` foreign of kind ``lazy``,
  ``pencil`` and ``posing`` of kind ``typing``; pure under the typing
  closure only. The test composes the question it answers,
  ``EigenPosing(OperatorPencil(A − S, F), K_MAP)``, and the inner inverse
  ``GreenOperator(A, S)``, and hands both to ``KEigenvalue``, so the foreign
  nodes are the claim's own inputs (the B3 specification, §5.2, read from
  the test's assertions);
* ``test_legendre_basis.py::test_gram_is_diagonal_and_exact_against_the_legendre_gauss_rule``:
  home ``quadrature``; ``basis`` and ``frame`` foreign of kind ``none``
  (``quadrature`` imports neither: the test builds the frame of a Legendre
  basis on a Gauss rule itself), its ``orpheus/geometry`` defs counted as
  an upstream layer;
* ``test_legendre_basis.py::test_endpoint_normalisation_is_exact``: home
  ``basis``, pure.

On that store the one support proposed runs against the order: the Gram
test over ``test_endpoint_normalisation_is_exact``, with ``c`` = 1.0
(``basis`` is not an ancestor of ``quadrature``).

.. _verification-cone-store-b4-placement:

The full-layer reading: B4's store
----------------------------------

When B3 landed no store of the whole numerics layer existed, and the
readings above, on samples of 10 and 180 tests, were not the layer's. The
store B4 captured at ``42b5ceca`` for its validation run
(:ref:`verification-cone-store-validation`) is the first. `[M]` 2026-11-07,
``scratch/verification_tree/cone_store/b4/runs/place.json``, its summary:

* 5 808 tests; 5 257 placed; 551 unplaced, of which 401 "no layer def"
  (B3's reason), 148 cone unknown and 2 not passed;
* pure under the three closures: 4 264 (runtime, the pure fraction 0.811),
  4 533 (lazy) and 4 941 (typing); the null's pure fraction 0.063, 0.060 and
  0.064 for seeds 1, 2 and 3;
* 993 straddlers and 20 ties; foreign nodes by kind, 1 119 ``lazy``, 1 736
  ``typing``, 320 ``none`` and 42 ``tie``; homes at 45 distinct labels;
* 798 declared ``rests_on`` rows, 138 at containment 1.0 and 513 below 0.8,
  and 275 edges unmeasured (221 "lower ran no def", 29 "upper not judged",
  25 "target not judged").

The census of the numerics layer reads this placement list by list (the
plan's ⏸ point of 2026-11-07). Its first list, the 401 tests with no layer
def, is what B5 answers (:ref:`verification-cone-store-reads`). The store was
captured before B5, so it holds no reads; the placement of the layer by
every channel is read on B5's re-capture
(:ref:`verification-cone-store-b5-acceptance`).

.. _verification-cone-store-place-gates:

The gates of B3
---------------

``tests/gates/tools/test_cone_store_phases.py`` gates the two cones on the
fixture tree of ``_cone_store_fixture.py`` with ``phases=True`` (a root
conftest that runs ``pkg.conf_dep.conftest_helper`` at import; ``test_p``
calls it again, ``test_b`` never does).
``tests/gates/tools/test_cone_store_placement.py`` gates ``place`` on the
fixture layer ``pkg/lay`` (``placement=True``), ranked by the fixture's
``PLACEMENT_ORDER`` so that rank order, alphabetical order and def count
disagree: ``beta`` imports ``zeta`` and ``s`` at top level, ``alpha``
imports ``beta``, ``gamma`` imports ``zeta`` lazily and ``beta`` under
``TYPE_CHECKING``, ``beta2`` imports nothing; ``pkg/up`` is the layer below
and ``pkg/top`` a layer above; ``tests/lay/conftest.py`` runs ``alpha.fd``
at import, so ``fd``, ``fb`` and ``fa`` are in every test's floor. Each
row's expectation is read off those sources by hand.

.. list-table::
   :header-rows: 1
   :widths: 10 58 32

   * - Row
     - Asserts
     - First red
   * - PH1
     - ``test_p``'s own cone holds ``conftest_helper`` although the
       grandparent ran it first; so does its selection cone
     - one recorder for both cones
   * - PH2
     - ``test_b``'s own cone lacks ``conftest_helper``; its selection cone
       holds it
     - the own cone written as the selection cone
   * - PH3
     - for every test, its own cone's nodes and the share's lie in its
       selection cone's; an edit to ``conftest_helper``'s body selects
       ``test_p`` and ``test_b``
     - ``select`` judging the own cones
   * - PH4
     - the grandparent's cone is stored as a set, holds
       ``conftest_helper``, and lies in every selection cone
     - the share not stored
   * - PH5
     - a store of format 2 is emptied when opened; ``own_cones`` and
       ``shares`` exist
     - ``FORMAT`` left at 2
   * - PL1
     - ``test_d_only`` (``alpha``, ``beta``, two defs of ``zeta``) is
       homed at ``alpha``; ``test_low`` at ``zeta``
     - a home by name or by def count
   * - PL2
     - ``test_d_only`` is pure; ``test_none`` adds ``beta2``, which nothing
       ``alpha`` imports: a straddler, ``beta2`` of kind ``none``
     - closure by rank order
   * - PL3
     - ``test_lazy_typ``: home ``gamma``, ``zeta`` of kind ``lazy``,
       ``beta`` of kind ``typing``; impure, impure, pure under the three
       closures
     - a walk that does not enter defs, or does not read
       ``TYPE_CHECKING``
   * - PL4
     - ``test_tie`` (``beta`` and ``beta2``, one tier): no home, the tie
       ``[beta, beta2]``, not pure
     - a tie broken silently
   * - PL5
     - ``test_low`` runs only ``zeta`` while the conftest ran ``alpha``:
       homed at ``zeta``, pure
     - placement read from the selection cone
   * - PL6
     - ``test_d_only``, whose every def is also in the floor, is placed at
       ``alpha``
     - placement by subtracting the floor
   * - PL7
     - ``test_layers``: home ``zeta``, pure, ``upstream_layers`` exactly
       ``{pkg/up/u}``, ``pkg/top`` under ``above_layer``
     - either read as a foreign node; the upstream count kept per layer
       root
   * - PL8 (3 ids)
     - a test that read nothing of the layer (B3: "no layer def"; since B5
       "no layer read"), a failed test and a test that starts a process
       are unplaced with their reasons
     - dropping them, or placing a partial cone
   * - PL9
     - ``c`` = 0.95 is proposed and 0.9 is not; an implied support is
       reduced away; a one-def cone is nobody's lower; cases are
       collapsed; equal cones are duplicates
     - ``>`` for ``≥``; no collapse; no reduction; no ``min_defs``; equal
       cones not grouped
   * - PL10
     - two runs give byte-identical JSON, with three null draws
     - an unseeded shuffle
   * - PL11
     - a store of format 2 is refused, before it is opened, by a message
       naming the own cones
     - a silent fallback
   * - PL13
     - ``ancestors`` on a 2-cycle, in three visiting orders, and on the
       live ``traced_memo`` pair
     - the memoised walk
   * - PL14
     - an import kind with no closure is refused by name
     - a kind read as ``none``
   * - PL15
     - each declared ``rests_on`` edge is measured on function cones:
       ``c`` = 0.9, 1.0 onto a case, 0 for a one-def lower, one edge for
       two cases, an unmatched target counted and not measured
     - no case collapse; ``>`` for ``≥`` in the counts; an unmatched target
       given a ``c``
   * - PL16
     - the null falls below the measured purity: `[M]` 2026-11-03 measured
       9/14, null 5/14, 4/14 and 7/14 (seeds 1, 2, 3), the draws distinct,
       on one denominator
     - the shuffle replaced by the identity (9/14 three times); one seed
       for every draw

The batteries (the test-architect, 2026-11-03; logs in
``scratch/verification_tree/cone_store/b3/battery_ph/arms.log`` and
``battery_pl/arms1.log``, ``arms2.log``): each of the five PH arms reddened
its target row (the share-count-only arm also PH3); on the placement rows,
the positive control (the selection cone read everywhere) reddened 8 of the
14 ids then present, and each of the 18 arms after it reddened its target
row. Two of those arms targeted PL12, the instrument row, retired with the
instrument region; PL13 to PL16 were written after that battery, each
with its first red in its docstring.


.. _verification-cone-store-reads:

Placing a test by what its run read (B5)
========================================

Increment B5 makes ``place`` home a test by what its own run READ, beside
the defs it ran (rulings 47 to 49 of the plan). Two pieces of code carry it.
``tools/cone_store/reads.py`` reads each executed first-party code object's
bytecode, resolves every load chain to the statement that binds its end, and
returns the reads of one test's run. ``tools/cone_store/placement.py`` sorts
those reads by channel and homes the test by a tiered rule. The store keeps
the reads in a new table, ``own_reads``, and its form moved to 4
(:ref:`verification-cone-store-store`). The specification is
``scratch/verification_tree/cone_store/b5/spec_b5.md`` (22 attacks, the
corrections C1 to C14, the gates and the acceptance); its prototype resolver
and probes are under ``b5/proto/`` and ``b5/probes/``. The reviews are
``b5/qa/qa_b5_review.md`` and ``b5/elegance/review.md``; the findings not
repaired in the build are #700.

Why the defs alone could not place 401 tests
--------------------------------------------

The placement of the B4 store (:ref:`verification-cone-store-b4-placement`)
left 401 of the layer's 5 808 tests unplaced with the B3 reason "no layer
def": their own cones held no def of ``orpheus/numerics``. The census's list 1
(``scratch/verification_tree/cone_store/b4/census/l1.md``, `[M]` 2026-11-07)
read each of them with two instruments: what the test ran (its own cone, by
pin kind and by package), and what it names (a static, module-qualified read
of the test function, its helpers and its fixtures, whose positive control
finds the placed home among the named nodes for 5 044 of 5 237 placed tests).
Each test falls in exactly one of four populations:

.. list-table:: The 401 tests of list 1
   :header-rows: 1
   :widths: 6 10 10 74

   * - Population
     - Tests
     - Functions
     - What the tests are, and the verdict
   * - G
     - 205
     - 47
     - Their own cones ran ``orpheus/geometry`` or ``orpheus/values`` defs
       only, and they name no numerics module. Their subjects (the subgroup
       lattice, G-spaces and quotients, the face-name bijection, geometry
       retirements, the memo key's encoding) moved below numerics in
       ``830b3713``, and the tests stayed. Misfiled; they are relocated per
       function when the folder is reorganised. B5 leaves them unplaced in
       this layer, and that is correct.
   * - S
     - 139
     - 68
     - Tests of a module's declared structure: a dataclass's fields, an
       enum's members, Protocol membership by ``isinstance``, a constant, a
       registry's population, a retired name's absence. They read a class
       statement or an assignment, which ran when the module was imported at
       collection, and they run none of the module's code. Correctly filed,
       and invisible to a placement by defs.
   * - A
     - 19
     - 19
     - Source audits: they open a numerics file as data (a ``DataPin``).
       Correctly filed, and invisible to a placement by defs.
   * - O
     - 38
     - 26
     - Tests of the test apparatus: oracle self-tests, fixture premises,
       the concept order's own reader, retirement gates. They name no
       ``orpheus`` module.

So 158 of the 401 (S and A) were blind spots of the instrument, not misfiled
tests. Population S is B4's second finding (a census over ``__subclasses__()``
that depends on which modules were imported, :ref:`verification-cone-store-validation`)
seen from the unplaced side. Why the own cone held nothing for them: the
memo's ``collect_pins`` mints a ``StatementsPin`` only for a file that has a
code object in the recording, and a test that reads a constant runs no code
of the constant's module. Ruling 47 (the user, 2026-11-07) answered list 1:
placement homes a test by what its run read, module-qualified, and reads a
source file opened as data as a node of the test. Ruling 48 answered
population O: a test of the test apparatus is a node below the gates that use
it.

The read relation, and why ``co_names`` is not it
-------------------------------------------------

A **read** of a test is a binding its own run loaded: a load chain of an
executed first-party code object (``m.a.b``: a load of ``m``, then the
attribute loads that follow it), resolved where the code's own namespace
binds the chain's root and then attribute by attribute, to the statement
that binds the chain's end. It is recorded as ``Read(path, name, channel)``:
the repository path of the module whose statement binds the name, the bound
name (``None`` for a read of the module object itself), and the channel the
read was reached by (below).

The design as first written (the plan's B5 design) took the names a code
object can read from ``_import_time.code_names``, that is ``co_names`` and the
identifier tokens of its string constants, the same names B2's entry-name
closure reads (:ref:`verification-cone-store-statements`). That cannot be the
read relation. ``co_names`` holds the globals a code object loads and the
attributes it loads in one unordered tuple, so it cannot pair an attribute
with its receiver: it cannot tell ``measure.WEIGHT_ATOL`` from
``quadrature.WEIGHT_ATOL``. Read literally ("an attribute name resolved in
every module the code reads"), it over-pairs. `[M]` the specification's
prototype on 658 tests of ``tests/gates/numerics`` (the 401 of list 1, 197
tests a proxy predicted to move, and 60 unmoved controls drawn with seed 5):
the design as written named a layer node that the corrected resolver does not
on 102 rows, and missed one it finds on 8. Correction C1 replaces it with
load chains read from the bytecode. For selection the over-pairing is safe,
since a larger cone only selects more; for a home it names the wrong node.

The bytecode reader
-------------------

``instructions(code)`` walks ``dis.get_instructions(code)`` once and returns a
``CodeReads``: the load chains of that one code object (never of the code
objects it holds), the names its import statements bind, and its star
imports.

* **A chain** starts at a load of a name: ``LOAD_GLOBAL``, ``LOAD_NAME`` or
  ``LOAD_FROM_DICT_OR_GLOBALS`` (scope ``global``), or ``LOAD_FAST``,
  ``LOAD_FAST_BORROW``, ``LOAD_FAST_CHECK``, ``LOAD_FAST_AND_CLEAR``,
  ``LOAD_DEREF`` or ``LOAD_FROM_DICT_OR_DEREF`` (scope ``local``). Each
  ``LOAD_ATTR`` that follows extends it. ``PUSH_NULL``, ``SWAP`` and ``COPY``
  leave it open; any other instruction closes it. ``CACHE``, ``NOP``,
  ``EXTENDED_ARG`` and ``NOT_TAKEN`` are skipped. Each chain keeps the line
  of its first load.
* **The import table** (C2). ``IMPORT_NAME`` takes its level and its
  ``fromlist`` from the two constants loaded just before it; each
  ``IMPORT_FROM`` adds a taken name; the store that follows binds a local
  (fast, cell or name) to an ``Imported(module, level, taken, fromlist)``. A
  load of that local resolves through the table before the module's globals.
  ``from m import *`` is kept apart as a star import. `[M]` the specification
  (``probes/inbody_imports.py``): 174 of 3 113 defs under
  ``tests/gates/numerics`` import ``orpheus`` inside their body, 153 of them
  test functions; without the table their reads name nothing.
* **A probe by constant** (C3). ``getattr(m, "X")`` and ``hasattr(m, "X")``
  load the callee, then the receiver, then the constant: when a global load of
  ``getattr`` or ``hasattr`` is followed by a chain and then a string
  constant, the constant extends that chain. Pytest's assertion rewriter
  stores such a constant apart from its receiver, so the string channel below
  also pairs the constants of the code with the modules it read.
* **A fused instruction is read as its two halves.** Python 3.14 joins two
  adjacent local accesses into one instruction (a superinstruction), whose
  argument is the pair of the two arguments. ``_operations`` yields each one
  as its two halves, and the table of them, ``_FUSED``, is derived from the
  opcode names rather than listed: an opcode whose name is two opcode names
  joined by an underscore. `[M]` 2026-11-07, Python 3.14.3: the derivation
  finds exactly four, ``LOAD_FAST_LOAD_FAST``,
  ``LOAD_FAST_BORROW_LOAD_FAST_BORROW``, ``STORE_FAST_LOAD_FAST`` and
  ``STORE_FAST_STORE_FAST``. An instruction the interpreter gains this way
  cannot fall through unread.

.. dropdown:: First got wrong: "the fused loads are two opcodes to list"
   :color: muted

   The first reader listed two fused opcodes and treated them as their own
   case. Two reviews found the two failures that shape has. qa's F2: the
   pair branch appended two bare chains and closed, so the ``LOAD_ATTR``
   after ``LOAD_FAST_BORROW_LOAD_FAST_BORROW`` (an attribute of the SECOND
   local, as in ``check(x, m.TOL)``) was dropped; `[M]` 9 such sites in 5
   numerics test files, none on a first-party module local. The elegance
   review's F8: ``STORE_FAST_LOAD_FAST`` was not listed, so
   ``import pkg as pk; value = pk.rd.high.h()`` on one line lost both the
   import's binding and the chain; `[M]` 0 sites in the tree, so the fixture's
   ``test_superinstruction`` is the only witness. Both are one defect: an
   enumerating guard over an open set. The repair splits every fused
   instruction into its halves before the walker sees it, from a table derived
   from ``dis.opmap``. Gates R18 (the store-and-load) and R21 (an attribute
   after a fused pair of loads).

**A cell is resolved through the def that binds it.** A name imported in an
outer def and read by an inner one is a cell. The resolution keeps, for each
executed code object, the import bindings of its cell variables, keyed by
(file, qualified name), and an inner code resolves a local through the tables
of the defs that enclose it (``f.<locals>.g`` is in ``f``) before its own.
The first build keyed the table by file alone, so the last executed code of a
file to bind a cell name won (qa's F3: ``helper_a``'s cell ``m`` read as
``helper_b``'s); gate R20 runs two such defs in one test.

**Which code objects are read.** Every code object of the own recording whose
file is first-party, except three kinds. A pytest hook implementation
(``pytest_*``) is the floor every test runs: `[M]` the specification,
``tests/conftest.py::pytest_runtest_makereport`` is in 5 808 of 5 808 numerics
own cones and reads a binding of ``orpheus/derivations``, so resolving it would
give every test an above-layer read (C4). A module body (``<module>``) that
ran during the test belongs to the imported module's provenance. A class body
that ran at import (a code object without ``CO_OPTIMIZED`` and with no
``<locals>`` in its qualified name) is the module's, not the test's.

Module-level and class-body code ran at import, outside the own recording, so
it is read another way: each first-party module the resolution touches is
compiled again from its bytes (``compile``, never ``ast.parse``), its module
body and class bodies are read with the same instruction reader, and each
chain is mapped to the statement that holds its line by the memo's
``_import_time.holding`` (the innermost statement whose lines hold it, so a
class's header holds its bases and decorators and its member statements hold
the body). The statements are the memo's own index of the file
(``traced_memo._source_index``), so there is one statement model, the
memo's.

Resolution to the statement that binds a name
---------------------------------------------

``Sources`` is the reusable half of a resolution: the interpreter's
first-party modules as they are after the test, indexed by real file, each
read once. ``resolve(sources, codes, test)`` is one resolution and returns a
``RunReads`` value (the reads, and per module path the names its pin is
projected on); it keeps no state between tests (the elegance review's F5).

``binding(module, name)`` answers where ``module.name`` is bound, in order:

1. A module dunder the import system binds and no statement binds
   (``__file__``, ``__name__``, ``__path__``, ``__spec__``, ``__loader__``,
   ``__package__``, ``__dict__``, ``__cached__``, ``__builtins__``,
   ``__doc__``, ``__all__``, ``__annotations__``) is a read of the module
   object (C8). So is a name whose value is a module.
2. **A re-export is followed by identity** (C7). When the import statement
   that bound the name (or, when no statement binds it, a star import) takes
   it from a module whose ``__dict__`` holds the very same object (``is``),
   the resolution continues in that module, and the module it passed through
   is recorded with channel ``reexport``: pinned, never a node. So a
   ``try: from x import y`` / ``except ImportError: y = …`` pair, and a
   ``TYPE_CHECKING`` import, stop at the module whose statement really binds
   the object. A re-export into another package's object is no first-party
   read. The reading code's own import of the chain's root is how the code
   reached the name, not a re-export, and is not recorded as one.
3. A non-import statement of the module binding the name ends the
   resolution there.
4. **A name the module holds that no statement binds** is an ``_Unbound``:
   a name set at run time (``setattr`` into a module, ``globals()[k] = v``).
   Its bare name joins the module's pin, so a statement that binds it later
   stales the pin; and, met on a node channel, it is recorded with channel
   ``unresolved`` and its dotted spelling (``pkg.rd.low.LATE``) as the name,
   counted, never dropped (C11). Met while following a binding's provenance
   (channel ``closure``), it is not the test's read, and only the pin keeps it
   (gates R19 and P9; the first build recorded it as the test's, the
   elegance review's F5).

A chain then steps attribute by attribute while the current end is a module;
an attribute a module does not bind is added to the tokens tried on that
module, and the chain stops.

The identity check is blind on an object Python shares: a small int, an
interned string, ``None``. `[M]` qa's probe: ``from other import N`` followed
by ``N = 3``, with ``other.N = 3``, resolves the read to ``other:N`` and keeps
``low:N`` only as a re-export. Gate R7b uses 99 and 5 and cannot see it; the
two candidate repairs are a design question (#700, §2).

The channels
------------

Each read carries a ``Channel`` (a ``StrEnum`` in ``reads.py``), whose
property ``is_node`` says whether the read makes a node of the test and
``is_pinned`` whether its name is one the module's pin is projected on.

.. list-table:: The channels of a read
   :header-rows: 1
   :widths: 16 56 14 14

   * - Channel
     - How the read was reached
     - A node?
     - Pinned?
   * - ``direct``
     - a chain of code the test ran
     - yes
     - yes
   * - ``decorator``
     - a chain in the decorators or defaults of a statement on the test's
       qualified name (its classes' headers, then its def), evaluated at
       collection, so the value reaches the test as an argument (C12)
     - yes
     - yes
   * - ``string``
     - an identifier token of a string constant (the docstring excluded) of
       code the test ran, held by a module that code reached as an object:
       the probe the assertion rewriter separates from its receiver (C3)
     - yes
     - yes
   * - ``module``
     - a chain that ends on a module object (``hasattr(m, "Old")``,
       ``m.__file__``), with no bound name; its pin is projected on the tokens
       the code tried on it, so binding a retired name later stales it (C8)
     - yes
     - yes
   * - ``via_apparatus``
     - a binding reached through the statement of a binding under ``tests/``
       (a case table, a fixture, a test-module constant), followed to the
       first binding outside the test tree (C6 as refined)
     - yes
     - yes
   * - ``closure``
     - a binding read by the statement of a binding outside the test tree
       (what built ``LEGENDRE``), followed one statement deep: the value's
       provenance
     - no
     - yes
   * - ``reexport``
     - a module a chain passed through by identity (C7)
     - no
     - yes
   * - ``ambiguous``
     - a binding whose file is in ``sys.modules`` under more than one module
       name (C10): counted beside the read, with no name
     - no
     - no
   * - ``unresolved``
     - a name bound at run time that no statement binds, spelled dotted (C11)
     - no
     - no

Ruling 48's apparatus and the refined C6 decide where the following stops. A
binding under ``tests/`` delivers a value, and its subject is the first binding
outside the test tree that its statement reads: ``CASES = (low.Rec(1.0),)``
read by a test reaches ``low:Rec`` as ``via_apparatus``, and the test is homed
at ``low`` (gate A2). A layer binding's own statement is provenance: what
built ``mid.BUILT`` is pinned as ``closure`` and is never a node (gates A3,
A4). `[M]` the specification's prototype: with no following at all, 65 of
the 139 S tests got a layer node; with following through apparatus only, 128,
of which 72 reached it through apparatus.

Not read, by declaration: a pytest hook implementation; module-level code
that ran during the test; a class attribute read through an instance whose
class the test never names (``box.LIMIT``, the instance from a fixture: gate
R15, one-sided); and a name built at run time (a string passed to
``importlib``: 2 S tests, `[M]` the prototype).

The statement pins of B5, and why they stale
--------------------------------------------

``cone_of`` resolves the reads in the worker after the test, for the own
recording only, under the recorder's ``unrecorded()`` (the resolution reads
sources and modules, and none of that is the test's). The reads go to the
store's ``own_reads`` rows. Beside them, one ``StatementsPin`` per module path
read enters the own cone, projected by the module's statement index on that
path's names: the names its pinned reads name (every channel but
``unresolved`` and ``ambiguous``), the tokens tried on it as an object, and
the bare names its provenance met unbound. The pin is the memo's own kind,
judged by the memo's projection rule
(:ref:`verification-cone-store-statements`), so it stales exactly when the
statements that bind those names change:

* an edit to a statement binding a read name (gate P2: ``TOL = 1e-12`` to
  ``1e-11``);
* a new statement binding a read name (P4: ``TOL = 2`` appended);
* a statement binding a name tried and absent (P6: ``GONE = 1`` added, the
  retirement gate's case);
* a deleted re-export (P5: the re-exporter's ``reexport`` pin);
* an edit to the builder of a value read through apparatus or as provenance
  (P7: ``Rec``'s field retyped);
* a statement binding a name the provenance met unbound (P9).

An edit to a binding the test did not read leaves its pin valid (P3, P3b):
that row tells a projection on the resolved names from a projection on the
run's whole entry-name closure, which would validate P1 too.

These pins are not the ones ``collect_pins`` mints. ``collect_pins`` projects a
``StatementsPin`` for every file with code in the recording on the run's
entry-name closure ``E′``, which is not what the test read, and keeps the names
in a set under the capturing checkout's ``.cache/sets``. So ``place`` never
reads names from pins: it reads the ``own_reads`` rows, which the store holds
(C9). One own cone may hold both pins for one file. No verb reads the B5 pins
today: they keep the own cone a freshness record of what the run read, and
the gates P1 to P9 are their consumers. Resolving adds nothing to the
selection cone (gate F3), so ``select``, which B4 validated, is unchanged. The
pin is minted in ``cone.py`` beside the memo's own mint site, because the
rulings froze the memo for B5 (#700, §4).

The tiered home
---------------

``placement.own_reads`` sorts a test's reads by region into an ``OwnReads``:
the layer nodes whose bindings it read on a node channel (``statement``, node
to names), the layer nodes whose files it opened as data (``data``, from its
``DataPin``\ s) and whether it opened a package file outside the layer too,
the apparatus (test-tree file to items), its reads of the package above the
layer, and its unresolved names. ``home_of`` returns how the test is homed as
one value of the sum ``Home``, which ``place`` matches once (the elegance
review's F1: a tag re-asked at five sites became a type):

.. list-table:: The tiers of a home (ruling 49), in the order they are tried
   :header-rows: 1
   :widths: 18 12 70

   * - ``Home`` member
     - ``placed_by``
     - When
   * - ``DefHome(nodes)``
     - ``def``
     - the test ran a layer def: the home is B3's, from the def nodes alone
   * - ``StatementHome(nodes)``
     - ``statement``
     - else, it read a layer binding on a node channel: the highest-ranked
       statement node, ties as for defs
   * - ``DataHome(nodes)``
     - ``data``
     - else, it opened layer files as data, not all of them: the
       highest-ranked data node
   * - ``ScopeHome("layer")``
     - ``layer``
     - else, it opened every node's file and nothing of the package outside
       the layer: no node is the home (C13)
   * - ``ScopeHome("package")``
     - ``package``
     - else, it opened every node's file and files of the package outside the
       layer too: a census finding, listed apart (C13, ruled)
   * - ``ApparatusHome(files)``
     - ``apparatus``
     - else, it read a test-tree binding or ran a test-tree def: its own file
       when it read there, else the one file, else a tie of the files (C14)

A test with none of these is unplaced with the reason "no layer read".

**Why tiered** (ruling 49, the user, 2026-11-07). Running a module's code is a
dependence on its behaviour; reading a value with none of its module's code
run is a dependence on an input. The home was designed for behaviour, so a
statement read never moves a def home: of the tests placed by defs, 0 change
home when the reads are added (gate V4 on the fixture; `[M]` qa's replay of
the B4 store through the build, 0 of 5 257 def homes and the pure fraction
0.811109 unchanged).

**The union home is kept as the refuted reading.** Every channel is reported
whatever homes the test: ``via`` (each layer node to its channels among
``def``, ``statement`` and ``data``), ``items`` (each node or apparatus file
to the names read there), and ``home_union``, the home every channel's nodes
would give together, highest rank winning. Union is refuted FOR the default
home; the fact it establishes is where a test's inputs come from. `[M]` the
specification's prototype: union moved 110 of the 197 proxy-predicted homes
(29 ``coefficient_layout`` to ``axis``, 22 ``measure`` to ``quadrature``, …)
and 0 of the 60 controls; ``test_exactness.py::test_claim_is_frozen``, for
one, builds an ``ExactnessClaim`` from ``LEGENDRE``, a constant of
``generating_measure``, and would move there. The summary reports
``homes_union``, ``union_moves`` (placed by def or statement, whose union
home differs) and ``statement_straddlers``: a def-placed test's statement
nodes its def home does not explain at runtime, each with the kind of the
first import edge that does (``statement_straddle``, read as ``foreign`` is).

**Purity and the null per channel.** Purity, the straddlers, the null and the
proposed ``executes`` supports are read over the def-placed tests, as in B3.
The statement-placed tests have their own ``pure_fraction_statement`` and a
null that shuffles the node label of each distinct (node, item) pair
(``null_pure_fraction_statement``); a module read is the item ``""``. A
fraction over no test is ``None``: the first build printed ``0.0``, which reads
as "every one straddles" where none was placed (qa's F6). `[M]` qa's capture
of 237 real tests: 12 statement-placed, 10 of them reading one node, measured
1.0 and null 1.0 on every seed, so whether the statement null separates on the
layer is an acceptance reading (AC7), not a fixture one (gate V6 records the
fixture's 23/24 against 24, 22 and 24 of 24).

Apparatus nodes and ``apparatus_supports``
------------------------------------------

Ruling 48 (the user, 2026-11-07): an oracle, a fixture or a case table is a
node of the DAG; its self-test homes there, and every gate that uses it rests
on it. C14 as ruled fixes the grain. The **apparatus node** is a test-tree
file (a helper module, a ``test_*.py`` that exports helpers, the test's own
module). Its **items** are bindings: the test-tree defs the test ran (its own
def, the code nested in it and the pytest hooks excepted, ``OwnDefs.helpers``)
and the test-tree bindings it read on a node channel. Why this grain: the
layer's nodes are files, so a home names a file either way; the premises
read module constants (``test_premise_the_metric_is_not_uniform`` reads only
``_G``, a module-level array of its own file), which a fixture-function grain
cannot express; and `[M]` the specification, 0 of the 153 test-tree defs other
test files run is a ``conftest.py`` fixture. A re-exporting test file is no
item (gate R16).

Each entry lists ``apparatus_items`` (``[path, name]`` pairs). The result's
``apparatus_supports`` proposes the premise edges: a test function ``P`` homed
at an apparatus file is the proposed lower of every other function ``U``
whose items at that file contain ``P``'s, with
``c = |I(P) ∩ I(U)| / |I(P)| ≥ threshold``, over functions (gate A6:
``test_uses_premise`` reads ``_G`` and runs ``low.f``, and rests on
``test_premise`` with ``c = 1.0``).

The declared ``rests_on`` edges: ``c`` and ``c_items``
------------------------------------------------------

The declared edges are measured over every test that passed with a known
cone, whatever homes it (``judged``), so that the population is the layer's
whole declared one. Each row now carries two containments: ``c``, over the
function's own first-party defs as in B3; and ``c_items``, over its defs and
its read items (a layer item spelled ``@<node>`` with its name, an apparatus
item by file and name). A lower function that ran no def and read items has
no ``c``: it is measured by ``c_items`` alone and counted under
``declared_items_only``. The cumulative ``declared_counts`` are over ``c``. An
edge whose lower ran no def and read nothing is counted under the reason
"lower ran no def and read nothing" (gate A9 measures an edge onto a premise
test by ``c_items``).

.. dropdown:: First got wrong: "the declared edges are measured over the placed tests"
   :color: muted

   The first build took ``judged`` to be the placed tests. In B3 it had also
   held the tests unplaced for "no layer def", so a lower test unplaced for
   "no layer read" (it ran only upstream or above-layer defs, or nothing)
   left the population, and its edge became neither a row nor a counted
   reason. `[M]` qa's replay of the B4 store (1 583 declared edges): the
   lower-level measurements fell from 1 019 to 902 while the reasons moved by
   6, so 117 vanished with no reason (qa's F1, an X2 failure: a filter
   upstream of the counter). ``judged`` is now every test that passed with a
   known cone. The B4 reading's 221 "lower ran no def" and the new key
   "lower ran no def and read nothing" are two predicates: the first counted
   every lower with no def, the second only those with no read either, the
   rest moving to ``declared_items_only``. AC10 re-counts the edges B4 could
   not measure, beside the defs-only count.

What was tried and refuted
--------------------------

* **Projecting a test's entry names onto the statements that bind them**
  (``l1_statements.py``, list 1's census) — refuted FOR homing. The entry
  names ``E′`` of a recording are not module-qualified (B2 accepts that a
  collision over-selects) and are one closure per recording, so the projection
  saturates: 45 to 47 candidate nodes on 80 tests, driven by ``__all__``,
  ``T`` and ``V``. The fact it establishes: every one of the 401 reads some
  numerics module-level binding by name.
* **The name tuple of a code object as the read set** (``co_names``): refuted for resolution (C1, above): it
  cannot pair a receiver with its attribute, and over-pairs on 102 of 658
  prototype rows.
* **Reading module-level and class-body code from the source** — ruled out by
  the recorder gate R7 (``test_cone_store_recorder.py``: no ``ast.parse``,
  ``ast.dump`` or digest of its own in ``tools/cone_store/``, one definition per
  quantity). An AST walk would be a second statement model beside the memo's
  index. The reader instead compiles the module's bytes, reads the module and
  class bodies with the same instruction reader, and maps each chain to its
  statement with the memo's ``holding`` over the memo's index. `[M]` the
  elegance review's p4: the span test the first build open-coded agrees with
  ``holding`` on 104 951 of 104 951 (statement, line) pairs over 229 files,
  and was retired for it.
* **The union home** — refuted FOR the default home (ruling 49, above); kept
  as ``home_union``.
* **No following at all** (the orchestrator's first wording of the closure
  ruling) — refuted FOR placing S:
  65 of 139 placed, where a test-module constant is how 72 S tests receive
  their subject.
* **Following provenance everywhere** — refuted for the home: the premise
  test would home at the builder of its constant (or at nothing, its builder
  being ``numpy``), and ``LEGENDRE``'s builder (``recurrence``, in
  ``orpheus/geometry``) would join every reader of ``LEGENDRE``.
* **Following a binding's provenance beyond one statement** — not built.
  `[R]` A layer binding's provenance is pinned, never a node, so its depth
  changes what the own cone's pins stale on and never a home. The deeper
  builders are import-time statements of modules the test's process imported,
  and the selection cone, which ``select`` judges, already holds them through
  the collection's share; B4 measured selection sound with it. One statement
  pins what the value read is spelled from. No measurement compares the two
  depths.
* **The fixture function as the apparatus grain** — refuted for this layer
  (C14, above).
* **Reading the read items from the statement pins**: refuted (C9): the
  pins of ``collect_pins`` name ``E′``, and their kept sets live outside the
  store.
* **The run-time module attribute as the resolver**: rejected as the primary
  route: a constant has no ``__module__``, and a class's can be reassigned.

Known limits
------------

The declared boundaries above (a hook, module-level code run during the test,
a class attribute read through an instance, a name built at run time) and the
review's findings filed rather than repaired, each measured at zero or one
site on the numerics layer today (#700):

* **Four reads dropped with no unresolved row** (qa's F4): a class-level
  case table read by a method's ``parametrize`` (``class T: CASES = …``; the
  class body's load resolves against the module's globals; `[M]` 1 site,
  ``test_registry.py::TestTheSelectorIsUNMOVEDByTheCarve``, def-placed, so its
  home is unaffected); a module-level ``__getattr__`` (0 in the tree); an
  in-function import of a module absent from ``sys.modules``; a file under
  two module names (0 ``ambiguous`` rows in qa's 237-test capture).
* **The identity check on a shared immutable** (above; #700, §2).
* **The test's own def is spelled two ways.** The worker hands the resolution
  the runtime identity ``(realpath(co_filename), __qualname__)``; placement
  re-derives the test's def from its id. The two differ on an inherited test
  method; `[M]` 0 of 1 382 test classes inherit one (#700, §3).
* **A second mint site for the statement pin** (#700, §4): ``StatementsPin``
  is minted in ``cone.py`` beside the memo's own site.

The gates of B5
---------------

``tests/gates/tools/test_cone_store_reads.py`` (groups R, P and F) and
``tests/gates/tools/test_cone_store_placement_reads.py`` (groups V, D and A)
are ``foundation`` rows on the fixture ``_READS`` of ``_cone_store_fixture.py``
(``build(root, reads=True)``): a layer ``pkg/rd`` ranked ``low`` below ``mid``
and ``rex`` below ``high`` (``READS_ORDER``), another package ``pkg/other``, the
apparatus ``tests/rd/_helper.py`` and a conftest with a hook and a fixture,
and one test per shape of read. Each row's docstring states its first red;
every expectation is read off the fixture's sources by hand.

.. list-table::
   :header-rows: 1
   :widths: 14 86

   * - Group
     - What it gates
   * - R (resolution)
     - each shape's exact ``(path, name, channel)`` rows: a module alias (R1),
       a dotted chain (R2), an in-body import (R3), a cell (R4, R20),
       re-exports through a module and through the package ``__init__`` (R5,
       R6), a fallback import and an import rebound by a later statement (R7,
       R7b), a name bound by a def and an assignment (R8), a dataclass and a
       Protocol (R9, R10), ``getattr``/``hasattr`` with a constant (R11), an
       absent name and a module dunder (R12, R13), a star import (R14), the
       instance-attribute boundary (R15), a re-exporting helper (R16), a
       ``parametrize`` decorator (R17), the fused instructions (R18, R21), and
       an unbound name met in provenance (R19)
   * - P (the pin)
     - P1 to P9 above: each B5 pin validates on the unedited tree and stales
       exactly on the edits listed; P8 round-trips ``own_reads`` and requires
       ``FORMAT`` 4 and a form-3 store emptied
   * - F (the floor)
     - F1: the conftest's hook reads ``pkg/other`` and no test but the one
       reading it itself holds that read; F2 its positive control; F3:
       resolving adds no ``DataPin`` or ``ListingPin`` to either cone
   * - V (placement by channel)
     - the tier (V1, V4), a re-exporter not a node (V2), a data home (V3),
       items (V5), the null per channel (V6), byte-identical JSON in-process
       and in two fresh interpreters under two hash seeds (V7), the union
       keys (V8), the above-layer statement count (V9), the unresolved count
       (V10), a form-3 store refused before it is opened (V11)
   * - D (the data home)
     - the ``layer`` home (D1), the ``package`` home (D2), data never above a
       def (D3)
   * - A (the apparatus)
     - a premise homed at its file (A1), a case table delivering its layer node
       (A2), resolution stopping at the first layer binding (A3, A4), a
       helper's reads being the test's (A5), the proposed premise edge (A6),
       a hook being no item (A7), the test not its own item (A8), a declared
       edge onto a premise measured by ``c_items`` (A9)

`[M]` 2026-11-07, the three files run serially under ``python -O -m pytest
-p no:cacheprovider`` on the uncommitted build: 74 passed, 0 failed, in
89 s; V11 landed after that run and passed in its file's re-run (23 of 23).
The new files hold 58 ids: R 22 (R1 to R21 and R7b), P 10 (P1 to P9 and
P3b) and F 3 in ``test_cone_store_reads.py`` (35); V 11, D 3 and A 9 in
``test_cone_store_placement_reads.py`` (23). The specification wrote 49; R7b,
R18 to R21, P3b, P9, F3 and V11 (a store of form 3 refused by a message
naming the own reads, before it is opened) were added after it, from the
battery and the reviews. ``test_cone_store_placement.py``
holds the other 17.

The B3 placement file is updated beside them: PL8's unplaced reason is now
"no layer read". The mutation battery of the specification (§3.7: each arm
a mutation applied to a copy of the built tree, the positive control C1's
literal pairing) ran on an earlier state of the build: `[M]` qa's review, 35
arms, 3 blind and each explained in its gate's docstring, and 3 of the then
52 rows (A5, P3 and V7) red under no arm. qa's three further arms reddened
A5, reddened P3b but not P3's first leg, and gave V7 a red that a broken
subprocess leg made unattributable. qa then found every battery reading void
on the tree that followed, where R18, V7 and V9 were red. The battery was
re-run on the final code in a separate worktree (`[M]` 2026-11-07, the
test-architect, ``scratch/verification_tree/cone_store/b5/battery/``): 41
arms, every one red; 58 of the 58 test functions of the two new files
reddened by at least one arm; the unmutated build 0 red; the positive
control (attribute chains dropped) reddens 23. Four arms of the earlier run
read blind and are not in the 41: F3 and F3b, replaced by F3c; C9b (the own
reads left out of the tables a store of another form drops) and C12b (the
first-party filter dropped from read resolution), which have no replacement,
so those two properties are not pinned by any gate (#700).

.. _verification-cone-store-b5-acceptance:

The acceptance on the numerics layer
------------------------------------

The specification's acceptance (§4, AC0 to AC10) ran on a re-capture of the
130 numerics test files at ``af5b78c8``, a store of form 4 holding 5 808
tests, against the B4 store as leg 1. `[M]` 2026-11-08: the capture took
7 h 28 min in 13 detached chunks
(``scratch/verification_tree/cone_store/b5/capture/``); the test-architect's
runner ``b5/acceptance/acceptance.py all`` checks that the placing code is
the capture's and writes ``acceptance.json`` and ``moves.md``. No check
failed.

* **AC0 and AC1, no def home moved.** 5 257 of the 5 257 tests B4 placed are
  in the new store and placed by defs, with 0 home drift against B4's
  placement; 0 of them has a home other than the def home recomputed from
  the store's own-cone ``DefPin`` rows.
* **AC2 and AC3, the unplaced populations.** 130 of the 139 S tests are
  placed by statement (the prototype predicted 128); the other 9 home at
  their own apparatus. The 19 A tests: 13 by statement, 2 by data, 2 at
  ``layer``, 2 at ``package``, as predicted. The 4 controls home at
  ``quadrature``, ``near_kernel``, ``quadrature`` and ``iteration``.
* **AC4, a record.** The home is among the nodes the static read names for
  142 of the 144 statement-placed tests; the 2 misses are the static
  instrument's (a re-export through ``spaces``, a helper reaching
  ``traced_memo``). Def-placed: 5 064 of 5 257. (B4's census counted 5 044
  of 5 237: its control left out the 20 tied tests.)
* **AC5.** 0 of the 205 G tests gains a numerics node (171 home at apparatus,
  34 are unplaced for "no layer read"). Of the 38 O tests, 36 home at
  apparatus, 1 by statement, and 1 is not homed
  (``test_the_invariance_module_is_deleted``, as predicted).
* **AC6, a record.** 116 tests have a union home other than their home, and
  229 def-placed tests are statement straddlers. Each move has a shape:
  a value delivered by apparatus, 54 (29 ``coefficient_layout`` to ``axis``
  through case tables of ``Axis`` values); an input value read from a higher
  node, 51 (21 ``measure`` to ``quadrature``, reading
  ``gauss_legendre_on_mu`` and ``NODE_ALIGNED``); a data read of a higher
  node, 6; a tie made by a shared read, 5. 0 of the 116 has no shape.
* **AC7, the null per channel.** The null is read over the tests a shuffle
  can move, those with two or more (node, item) labels: a test with one label
  stays pure under every shuffle, so over all tests the null of a channel of
  mostly single-label tests reads near its measured value and says nothing.
  The def channel, 5 026 movable of 5 257: 0.802 against 0.020, 0.017 and
  0.021 (seeds 1, 2, 3). The statement channel, 71 movable of 144 (80 labels
  over 28 nodes): 0.972 against 0.141, 0.127 and 0.056. Over all tests the
  readings are 0.811 against 0.060 to 0.064, and 0.986 against 0.535 to
  0.576. The data channel holds 2 tests and has no null. The runner's null
  equals ``place``'s on every seed of both judged channels.
* **AC8 and AC9.** 0 unresolved reads, 0 ambiguous files; 0 tests hold the
  ``pytest_runtest_makereport`` read (the floor).
* **AC10, a record.** Of 976 declared ``rests_on`` rows, 798 are measured
  over defs and 178 only through read items (85 of them at 1.0); the reason
  "lower ran no def" (221 in B4) is now "lower ran no def and read nothing",
  43.

The placement of the layer by every channel: 5 257 by defs, 144 by
statement, 216 at apparatus, 2 by data, 2 at ``layer``, 2 at ``package``
(5 623); unplaced 148 cone unknown, 35 no layer read, 2 not passed (185).

The supports against the order on this store, read for the census (its list
2, ``scratch/verification_tree/cone_store/b5/census/l2.md``), are what led
to ruling 50 (the subject's containment, above). `[M]` 2026-11-08, ``place``
re-run on the same store with every other key of the report unchanged.
Of the 6 784 supports along the order, the subject's containment fails for
1 161 (for 58 of them the upper ran none of the subject); 1 805 are dropped
in all (the 1 161, and 644 now implied through another support) and 1 605
added, 6 584 in all. Against the order, 126 to 79: 45 kept, 81 dropped, and
34 that the transitive reduction had hidden behind middles the ruling
removed. ``support_counts`` at 0.95 fall from 47 312 to 39 598, and the
near-duplicate groups from 105 to 76. In 79 of the 79 remaining pairs the
upper test's own source names the lower's home node (the same check reads
False for 57 of the 81 dropped pairs): the test runs a node its home does
not import, a straddle. That is a reading of the tests; it does not exclude
an edge from the upper's home module that the import graph misses (a
dynamic import, a dispatch).


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
  ``filterwarnings`` or ``addopts`` setting (qa's review, U4; gate S14);
* the bytes of the capture's instrument (``environment.INSTRUMENT``): the
  memo's own sources (its ``_MEMO_SOURCES``: ``traced_memo.py``,
  ``_traced_memo_boot.py``, ``_import_time.py`` and
  ``orpheus/values/content.py``) and every module of ``tools/cone_store``
  except the two that only read a store, ``selection.py`` and
  ``placement.py`` (a new module counts until it is named a reader). A cone
  holds no pin on the instrument (the memo-source pin is dropped, and the
  apparatus is never recorded), so a store is valid only under the
  instrument that captured it: the memo's own rule for its entries, read as
  a key. An edit to the instrument selects every test. That is a cost: an
  edit to ``store.py``, ``plan.py`` or ``__init__.py`` of
  ``tools/cone_store`` invalidates the whole store, and every test is
  captured again.

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
   kind (T1; a ``StatementsPin`` is contained when the store observes every
   statement of its file that the alone run observes, never by its names, as
   the Gotchas below explain); the same on a real-suite sample chosen by mechanism, 14 tests
   in five files plus one spawner asserted unknown (T2); the same on sample
   N's 200 tests (T3, ``slow``); the capture's and the alone run's outcomes
   agree (T4). The alone capture shares the recorder with the store, so the
   comparison is independent on the route axis only (X4); the recorder's own
   ground is R1 to R5, whose expected sets are written by hand from the
   fixture's source.

The gates, ``tests/gates/tools/test_cone_store_*.py``, every one
``foundation`` (`[M]` 2026-11-04, ``python -O -m pytest --collect-only`` over
the 10 files: 129 ids, and 0 of them left by ``-m "not foundation"``; 4 of them carry
``slow``, M1, T3, SG4 and SG9, none parametrised; and the recorder's
apparatus gate beside them, ``tests/gates/numerics/test_traced_memo_apparatus.py``,
21 ids by the same command):

.. list-table::
   :header-rows: 1
   :widths: 26 14 60

   * - File
     - Rows
     - Claims
   * - ``test_cone_store_recorder.py``
     - R1 to R9
     - exact pin sets against hand-written expectations; every digest is
       the memo's; the store defines no digest and parses no source (R7),
       installs no ``sys.monitoring`` tool or audit hook of its own, both
       of a child's recordings being the memo's ``Recording`` (R8), and
       asks ``fork_blockers`` rather than re-spelling the guard (R9)
   * - ``test_cone_store_collection.py``
     - C1 to C6, H1
     - a collection fill is attributed to its module; another module's fill
       cannot hide it (C2); the grandparent's code reaches every test (C4);
       a package ``__init__`` reaches every test under it (C5); the rows do
       not depend on module order (C6); the in-suite monitor's blindness to
       a fill is asserted, not assumed (H1)
   * - ``test_cone_store_selection.py``
     - S1 to S24, P1 to P3
     - ``select``'s rule, the three controls, the environment key per
       component (S5, S6), the shadow pins (S12), a new conftest selects
       every test (S13, S19), the ini bytes (S14), a collection error (S15),
       a read-then-rewritten file (S16), an edit during the capture (S17),
       an off-path conftest's session hook (S18), a user
       ``pytest_ignore_collect`` stays recorded (S20), a ``pytest11`` entry
       point under an unchanged version (S21), ``PYTEST_PLUGINS`` (S22), a
       reader of its own memo store selected by an edit to the generated
       code and by nothing else, its own cone without the generation's
       defs (S23, ERR-125), and a def such a reader also ran itself kept in
       its own cone (S24)
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
     - SG0 to SG10
     - B2's statement grain at the store: no edit selects nothing (SG0);
       control 1 exact on a fixture (SG1, SG2) and on a clone of the
       production tree (SG4, ``slow``); recall at the store for an effect,
       a re-export, a default argument and qa's S2a and S2b″ (SG3); a store
       in B1's form is emptied and every test is new (SG5); an import that
       starts raising stays selected after a ``select`` (SG6); a module
       under another ``sys.path`` entry is judged in the store's spelling
       (SG7); the capture's baseline covers every repository source (SG8);
       the capture pins none of its own code, in either cone, on the
       production tree (SG9, ``slow``); a tree package shadowed by the repository's own is
       measured in the run's precedence (SG10); a guarded probe of a name
       new to its module selects the prober, with a precision leg (U2)
   * - ``test_cone_store_phases.py``
     - PH1 to PH5
     - B3's two cones and the grandparent's share
       (:ref:`verification-cone-store-place-gates`)
   * - ``test_cone_store_placement.py``
     - PL1 to PL11, PL13 to PL16
     - ``place`` on a fixture layer
       (:ref:`verification-cone-store-place-gates`)
   * - ``tests/gates/numerics/test_traced_memo_apparatus.py``
     - 6 tests, one of them 16 arm cases (21 ids)
     - the recorder never records its apparatus, on every hook path and
       in a nested recording, and records the rest of its file
       (:ref:`the recorder's apparatus <verification-cone-store-apparatus>`)
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


.. _verification-cone-store-validation:

The validation run (B4)
=======================

**The question.** Does ``select`` list every test an edit can redden? A
selection that omits such a test is an under-selection, the one failure the
store must not have: the test would be skipped while it fails. An
over-selection only costs time. The plan's ruling 41 asks the question of
sampled mutants, each run both in full and through the store; ruling 45 (the
user's, 2026-11-04) fixes the terms below. The spec, the driver and every
per-mutant record are under ``scratch/verification_tree/cone_store/b4/``
(``spec_b4.md``, ``driver/``, ``runs/``).

**The terms.**

* **The tree and the store.** The run is at ``42b5ceca``, against a store of
  the numerics layer captured at that commit in a separate worktree.
  ``select`` lists every test stored under another environment key, so the
  store and the tree it judges must be one commit.
* **The ground truth is the full-suite context.** Every full run collects
  the whole ``testpaths`` tree, as the merge gate does, so all 12 conftests
  and every test module are imported, and only the numerics tests run
  (``--deselect`` of every other path; a control checks 5 808 kept, equal to
  the layer collected alone, and 12 conftests loaded). The store promises
  the verdict a full run gives (ruling 42), and a run of the numerics folder
  alone loads one conftest.
* **A mutant is an edit on disk** (``select`` judges file bytes): an exact
  replace asserted to occur once, the pristine file copied aside first, and
  each restore witnessed three ways: the bytes, ``git status`` against its
  reference, and a pristine ``select`` equal to the no-edit selection. Edits
  to the memo's own source and to ``tools/cone_store/`` are excluded: they
  move the environment key, which is another claim.
* **A red** is a test whose outcome differs from the pristine baseline. A
  module that fails to collect counts its tests as red.
* **A counted mutant** is one that can miss: ``select`` exits 0, lists fewer
  than every test, gives no environment-key reason, and the full run reddens
  a test outside the no-edit selection that is neither new nor in a module
  that failed to collect. The acceptance asks 0 misses over at least 60
  counted mutants that kill, which bounds the per-mutant miss rate below 5%
  at 95% confidence.
* **The frame** is 85 mutants: 49 def bodies spread over the concept
  order's tiers (stage A, whose miss would stop the batch), and rows by
  shape: import-time constants, default arguments, dataclass fields,
  ``values/`` and ``geometry/`` code, test helpers, ``parametrize``
  constants, a conftest, a new file, a collection error, nulls (a docstring,
  an unread binding, an assert stripped by ``-O``) and controls.

**The result** (`[M]` 2026-11-05 to 2026-11-06, serial, ``-O``, 42 h 14 min):

* **0 under-selections over 61 counted mutants that killed**, so the
  per-mutant miss rate is below 3/61, about 4.9%, at 95%. A positive control
  of the miss counter (one red id dropped from a copy of a selection) read
  its miss. 24 mutants did not count: 13 whose reds all lay in the no-edit
  selection, 4 that cannot miss (``select`` refused, or listed every test),
  and the 2 floor rows, 3 nulls and 2 controls.
* **The controls.** The positive control (the root conftest's report hook
  marking every pass failed) selected exactly the layer's 5 808 tests and
  reddened more than 5 000. A syntax error in a module the root conftest
  imports made ``select`` exit non-zero instead of listing a short
  selection. Three pristine full runs agree, and no test was flaky.
* **The restores**: 85 of 85 witnessed on all three legs.
* **Precision.** A mutant that killed selected a median of 72 tests beyond
  the 150 selected with no edit (the cone-unknown and the not passed), out of
  5 658. The ratio of tests selected beyond the no-edit selection to tests
  reddened beyond it has a median of 8.3, from 1.0 to 592.
* **The selection run alone.** Each mutant's selection was also run as a
  user runs it (``pytest <ids>``), and its verdicts agree with the full
  run's in 77 of 77 mutants. The one flag the record raised, on the mutant
  ``T07b``, is an artefact of the record: both runs hold the same collection
  error of ``test_inverse_metric_operator.py``, keyed by the file's path,
  which the full side's filter on selected ids drops.
* **Placement predicts propagation.** Over the 59 mutants analysed, no red
  test contradicts its placement: each is homed at the mutated node or
  downstream of it, is a straddler whose foreign nodes include it, or is
  unplaced. One red test's own cone lacks the mutated node altogether, a
  dependence of another kind (the second finding below).
* **The cost.** A mutant's ``select`` took a median of 248 s and its full
  run 874 s (#699).

**Two findings.**

* **A null control predicted wrong.** A docstring edit in
  ``orpheus/numerics/generating_measure.py`` was predicted to select only the
  no-edit selection. It selected 8 tests more and reddened none, and all 8
  carry the reason "data … generating_measure.py changed": they read the
  module's source as data (the concept-order census, import and AST audits),
  so their verdict can depend on its text, and selecting them is correct.
  The corrected predicate, the no-edit selection plus the tests holding a
  ``DataPin`` on the file, holds for 8 of 8.
* **A dependence through import, not through a def.** Under a syntax error
  in ``orpheus/numerics/near_kernel.py``,
  ``test_axis_identity.py::test_K1_equality_and_hash_are_defined_once_on_the_base``
  went red, and its own cone holds no ``near_kernel`` def. The test is a
  census over ``__subclasses__()``, which sees the classes of the modules
  that were imported; importing a module runs its module code, which the
  cone records as statements, never as defs. ``select`` covers the test
  through the module's import-time statements. Placement by defs cannot see
  such a dependence; since B5 a test that reads a module's bindings without
  running its code is homed by those reads
  (:ref:`verification-cone-store-reads`).

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
  names whole (gate row U2, now a recall row; at the store,
  ``test_cone_store_statement_grain.py::test_u2_a_guarded_probe_of_a_name_new_to_its_module_selects_the_prober``,
  with a precision leg: a new statement binding a name the run never names
  selects nothing).
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

One further limit is a property of the evidence, not of the store:

* the capture's cost on the whole suite has not been measured. Per child,
  ``collect_pins`` takes 1.8 to 1.9 s, mostly re-reading the inherited
  floor, and the own cone 0.10 s more (`[M]` 2026-11-03,
  :ref:`verification-cone-store-own-cone`); `[R]` the numerics layer's
  5 780 tests take about 2.6 h. The cost is paid once per environment and
  then only for stale tests; amortising the floor's read per collecting
  parent is #696.


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
* **The own cone is not the selection cone minus the floor.** A def the
  conftests ran is in a test's own cone when the test runs it too. A
  reading of "what this test exercises" computed as the selection cone less
  the grandparent's share drops every such def; read ``own_cones``
  (``ConeStore.own_pins``).
* **Compare the two cones by node, never by row.** One child's two cones
  name a module by different ``StatementsPin`` and ``ImportOutcomePin``
  rows (each is the projection of its own run), so a row-level containment
  test between them fails on a correct store.
* **Compare two statement pins of one file by the statements they observe,
  never by their names.** The entry names a ``StatementsPin`` keeps are one
  set per run, over every file the run's code came from, and the pin's digest
  projects that set onto one file's import-time statements. Two runs of one
  test can therefore keep different name sets while observing the same
  statements of a given file: the set depends on what else ran in the
  process, so on the test order (B3's acceptance, plan row X1.1, saw it gain
  the name ``fuzzy`` with the order). A name in one set and not the other
  that binds nothing in the file changes no statement either pin observes.
  A containment check between two such pins compares the observed
  statements, ``_import_time.selection`` over the file's statement index
  with each pin's kept names, as ``_misses`` in
  ``tests/gates/tools/test_cone_store_ground_truth.py`` does; comparing the
  name sets reads a miss where nothing is under-selected. `[M]` 2026-11-03
  (``f78085d8``): two T3 runs failed comparing names, one of them on the
  tree before the over-selection fix (``9195dcf5``); in 10 of 10 printed
  pins the alone run's projection equalled the store's digest for the same
  test from a passing run. On the kept stores of a passing T3 run the
  statement predicate read 0 misses over 200 tests; 44 000 alone pins
  observe some statement by name, beyond those every run observes, so the
  predicate can fire (T3 asserts this activation leg); and with the store's
  names emptied as a control, 190 of the 200 tests read misses.
* **Opening an old store empties it.** ``ConeStore`` and ``open_store`` drop
  the tables of a store of another format when they open it; only
  ``place`` reads the format first (``stored_format``). Copy an old store
  before opening it with this code if its rows are wanted.
* **A home is the most downstream layer code the test ran, which may be an
  input of the claim rather than its subject.**
  ``test_keigenvalue_recovers_dominant_eigenvalue`` asserts ``KEigenvalue``'s
  dominant ``k`` (``iteration``) and is homed at ``green_operator``, because
  it builds the Green operator it hands in. Read a home as "the most
  composed code this test runs", and a straddle of kind ``none`` as a test
  that composes two concepts itself.
* **The audit hook runs untraced.** CPython calls an audit hook with
  ``sys.monitoring`` paused unless the hook sets ``__cantrace__``, so code
  that runs only inside ``_audit`` is never started under any recording. A
  gate that wants to see such code recorded must set ``__cantrace__`` on the
  hook, as the apparatus gate's traceable run and arm cases do, each in a
  fresh interpreter so the traceable hook never reaches a recording around
  the test.
* **Only a process's first recording installs the probes.**
  ``Recording.start`` replaces ``os.stat`` and ``os.lstat`` with probes
  only while no recording is active, so a gate that needs ``_probing`` to
  run must start the first recording of its process. Under the cone store's
  capture every test already runs inside a recording; such a gate runs its
  check in a fresh interpreter (`[M]` 2026-11-03: the apparatus gate's
  ``[_apparatus]`` and ``[_probing]`` arm cases, green plainly and red
  under capture until they moved).
* **The summary's placed count is of the tests placed by their defs.** Since
  B5 a test is also placed by its reads, its data or its apparatus, and
  ``placed_by`` counts each tier. Purity, the straddlers, the null and the
  proposed ``executes`` supports are read over the def-placed tests; the
  statement-placed tests have their own pure fraction and null. Reading
  ``placed`` as "every test with a home" undercounts.
* **A home read from a value is not a home read from code.** A test homed by
  ``statement`` read a binding of that node and ran none of its code: its
  claim rests on the node's declared structure or on a value the node
  delivers, never on code of the node that its run executed. Read
  ``placed_by`` before reading a home.
* **The concept order is loaded from a test module by path.** ``--order``
  names a file, which ``__main__`` loads with
  ``importlib.util.spec_from_file_location``; that works because the
  order's module uses absolute imports only, and a relative import added
  to it would fail there and not in pytest.


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
   * - 2026-11-08
     - **Rulings 51 to 53: a declared ``rests_on`` is an execution
       support**, held by the support predicate (``c`` and ``c_home``), and a
       class or file target holds for an upper that contains some lower
       under it. Census list 3: 170 of 443 declarations on the numerics
       layer hold, 273 fail (the marker repair is owed).
     - #358, #701
     - ``docs(verification): census list 3, rulings 51 and 52``;
       ``fix(cone_store): place reports each declared edge's verdict``
   * - 2026-11-08
     - **Ruling 50: a support contains the lower's subject.** The census's
       list 2 found 62 of the 126 supports against the order whose upper ran
       none of the lower's home-node defs. ``_supports`` now also requires
       ``c_home ≥ threshold`` over the function's own defs at its home nodes
       (``OwnDefs.at``); on the numerics store this removes 1 161 of the
       6 784 supports along the order. Support records and declared rows
       carry ``c_home``; one ``_contained`` measures every containment.
     - #358
     - ``fix(cone_store): a support contains the lower's subject``
   * - 2026-11-07
     - **Increment B5: placement by what a test's run read** (rulings 47 to
       49). ``tools/cone_store/reads.py`` reads each executed first-party
       code object's load chains from its bytecode (fused instructions split
       by a table derived from ``dis.opmap``), resolves each to the statement
       that binds its end (re-exports by identity, unbound names counted),
       and sorts the reads into channels; the store's ``own_reads`` table and
       ``FORMAT = 4``; one ``StatementsPin`` per module read in the own cone;
       the tiered home as the ``Home`` sum (def, statement, data, layer or
       package, apparatus) with the union home reported and refuted for the
       default; apparatus nodes and ``apparatus_supports``; ``c_items`` on the
       declared edges, over the restored ``judged`` population (qa's F1).
       The motive is list 1's 401 unplaced tests. Review findings not
       repaired: #700. The acceptance on a re-capture of the numerics
       layer (AC0 to AC10) passed, 2026-11-08.
     - #358
     - ``af5b78c8``
   * - 2026-11-03
     - **B3's acceptance: a no-edit select over-selected** (ERR-125). Plan
       row X1.3 counted 229 numerics tests where 149 were owed. 70 readers
       of a memo store under their own ``tmp_path`` held a ``ChildPin``
       judged against the default store; a child of another store is now
       pinned by the pins its lookup validated it on (``Hit.pins``,
       ``ReadChild``), each source path placed where the read found it,
       and the limit "a test reading children from its own temporary memo
       store is always selected" is retired. 10 source-walking tests held
       a ``__pycache__`` listing, which is no longer observed. Gate S23,
       and T2's reads leg. T3 then compares statement pins by the
       statements they observe (the Gotchas).
     - #358
     - ``b65674fc``, ``f78085d8``
   * - 2026-11-06
     - **The validation run (B4) on the numerics layer: 0 under-selections
       over 61 counted mutants that killed** (the bound below 5% at 95%), the
       selection run alone agreeing with the full run in 77 of 77 mutants;
       the merge gate stays the full suite (ruling 46).
     - #358
     - ``42b5ceca`` (measured), plan ``2c2cea2c``
   * - 2026-11-04
     - **The numerics re-capture: X1.3 met, and an own cone that held
       another process's defs.** A no-edit ``select`` at ``f78085d8`` listed
       150 of the layer's 5 808 tests, its 148 cone-unknown and 2 not passed,
       with nothing stale. Placement then moved 70 tests, exactly the
       readers of a store under their own ``tmp_path``, from ``traced_memo``
       to ``convergence``: the pins of a foreign child had entered the own
       cone. They are now ``CollectedPins.through_reads``, in the selection
       cone and out of the own cone unless the run ran them too
       (``ReadCones``); S23's own leg, and S24 for a def both ran.
     - #358
     - ``1e989da7``
   * - 2026-11-03
     - **Increment B3: two cones per test, and placement.** The own cone (a
       second ``Recording`` around the test's protocol, table
       ``own_cones``) beside the selection cone (``cones``, unchanged); the
       grandparent's cone stored once in ``shares``, from which ``status``
       derives its count (the ``grandparent_pins`` meta is gone);
       ``FORMAT = 3`` and ``stored_format``; the worker strips its own code
       from both cones. ``place`` (``tools/cone_store/placement.py``): the
       home, the tie and purity under the three closures of the concept
       order's import graph (``import_kinds``, the gate's walk, now with
       kinds and an exact ``TYPE_CHECKING`` predicate), the ``Region`` sum,
       the unplaced reasons, the seeded null, fractional-containment
       supports reduced over the condensation; ``ancestors`` by
       reachability (the memoised walk failed on a cycle). Gates PH1 to
       PH5 and PL1 to PL11, PL13 to PL16; declared ``rests_on`` edges
       measured on the function cones. The recorder's apparatus is declared
       and never recorded (``_apparatus``), after the archivist found its
       hooks in every own cone and the instrument fallback homing 6 of 10
       sample tests at the recorder; the fallback and PL12 were retired,
       and the environment key digests the capture's instrument.
       Follow-ups: #695, #696.
     - #358
     - ``9195dcf5``
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
     - ``6e52b3cf``
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
