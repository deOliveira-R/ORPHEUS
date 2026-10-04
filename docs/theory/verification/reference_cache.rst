.. _verification-reference-cache:

=============================================================================
The reference cache: a traced memo keyed on what ran
=============================================================================

.. Machine header. Owner of: the traced memo (#405 P3), the primitive
   that memoises a reference reading on disk under a key that misses
   exactly when something that produced it has changed; why the identity
   of a generator is its recorded execution and not a declared version or
   a static call graph; why a miss is generated in a fresh process; the
   exact key; the manifest and its eleven pin kinds; validation that runs
   nothing; the exact payload; the process boundary; the bypass; how the
   memo relates to the reference certificate; the exact medium's
   exclusion; the clients; the limits no manifest can see; the review
   round's defects and their fixes; the ambient-state step (the slab
   route as an argument, the mpmath precision as a local context).
   Not owned here: the ontology of a reading, a claim and a certificate
   (readings_and_certificates.rst, verification-reference-architecture);
   the content encoder's canonical forms
   (structured-geometry-content-identity); the trajectory resolvent's
   reading (references/trajectory_resolvent.rst,
   trajectory-resolvent-reference-reading); the slab routing's physics
   (references/peierls_nystrom.rst, theory-peierls-slab-polar-g5-routing).
   Plan of record: .claude/plans/reference_cache.md, from "P3 opened
   (2026-10-04)" to its end, and discussion 2; specification of the
   gates: .claude/plans/reference_p3_spec.md. Written at P3's close,
   2026-10-04, against branch refactor/p3-ambient-state at 39ee20f2;
   updated for the second review round at 994ba740 and for the generator
   contract and the third review at bacf0787.

A reference reading costs seconds to minutes. The multi-region sphere
of the A|B|A cross-check, at the coarse resolution the gates use, takes
about four seconds to solve; the trajectory-resolvent reference rows,
the S\ :sub:`N` cross-checks that read it and the multi-region solver
tests, the fourteen files of the specification's baseline, summed to
5 327 s, about 89 minutes, per pass before this phase (`[M]` 2026-10-04,
one run per file, concurrent with other work, so an upper bound; seven of
the files are in the table of "The numbers" below). A reading changes
only when something that produces it changes, so the same answer was
recomputed on every run of the suite. This chapter documents the
primitive that stops that: :class:`~orpheus.numerics.traced_memo.TracedMemo`,
a pure function of content-identified arguments, memoised on disk under
``.cache/references/`` (ignored by git), whose entry is served again for
as long as nothing the call depended on has changed. Its dependencies are
not declared: they are **recorded**. On a miss the call runs in a fresh
interpreter that records every function that starts, every file it opens
or probes and every directory it lists, and the entry keeps that record, the
**manifest**, beside the answer. A later lookup compares the manifest with
the checkout, without importing or running anything, and serves the answer
only if every recorded dependency is unchanged.

The memo raises the development cadence; it certifies nothing. A wrong
reference served from the cache is exactly as wrong as the same reference
computed, and the reference certificate of
:ref:`verification-reference-architecture` sits above the memo, untouched
by it (the section "The memo and the certificate" below).

.. contents::
   :local:
   :depth: 2


Key facts
=========

* **The identity of a generator is its recorded execution**: the
  normalised source of every first-party function that ran, the skeleton
  of every first-party module whose code ran, the version of every
  third-party distribution whose code ran or whose file was read, the
  interpreter, every file read or probed and every directory listed
  (outside the import system), the declared environment, the starting
  working directory when a relative path was read, every memo entry read
  (checked recursively), and the memo's own source (the user's ruling of
  2026-09-24, discussion 2 of the plan). A declared version and a static
  call graph were refuted for this question (the section "Why a recorded
  execution" below).
* **A miss is generated in a fresh process** (the user's ruling of
  2026-10-04, architecture B). A process born for one call has nothing in
  memory, so a ``functools.cache`` hit, a ``cached_property`` computed
  earlier and a monkeypatch cannot hide a dependency from the record.
  In-process memoisation in ``orpheus/`` stays as it is. The price is one
  interpreter start per miss, about 0.65 s (`[M]` 2026-10-04, three runs).
* **The key is exact, finer than content equality.** It is
  :func:`~orpheus.numerics.content.encode_exact` of the function's
  identity, its signature-bound arguments (defaults applied, each declared
  canonical form put in) and the platform tag. ``8`` and ``8.0``, ``-0.0``
  and ``0.0``, an integer array and its float twin are different keys,
  because a function can tell them apart even though ``==`` cannot; and
  the key is exact all the way up: an object is its constructor form, a
  container carries its concrete type, a mapping its order.
* **Validation runs nothing.** A lookup re-hashes the manifest's pins
  against the current files and returns one of four verdicts,
  :class:`~orpheus.numerics.traced_memo.Hit`,
  :class:`~orpheus.numerics.traced_memo.Absent`,
  :class:`~orpheus.numerics.traced_memo.Stale` (with every reason) or
  :class:`~orpheus.numerics.traced_memo.Corrupt`. Only a ``Hit`` is
  served; every other verdict regenerates.
* **The payload is exact and is never a pickle**: floats as
  ``float.hex``, arrays in the entry's ``.npz`` with their dtype and
  bytes, numpy scalars with their dtype, tuples, and only the dataclasses
  the function's return annotation names, rebuilt through their
  constructors. A subclass of an admitted type is refused, never written
  as its base. An entry is one file, replaced atomically.
* **The memo sits below the certificate.** It caches what a derivation's
  ``evaluate`` returns, an establishment or an uncertified value, never a
  certificate or a verdict; ``ReferenceSolution`` re-checks every claim
  against that evaluation, cached or computed.
* **Three clients in P3**: the two multi-region trajectory-resolvent
  solvers (``solve_greens_function_sphere_mr`` and
  ``solve_greens_function_cylinder_mr``), and the trajectory-resolvent
  derivation's ``evaluate``, a memoised method keyed on the derivation's
  content. The exact infinite medium is **not** a client: its reading
  costs 0.01 s, less than one interpreter start (the user's ruling, Q5).
* **A test that monkeypatches anything a generation runs reads under**
  ``with bypass():``, because a patch never reaches a fresh process.
* **A generation may start only a declared machine query** (``uname``)
  and sees only a declared environment; every other started process is
  refused.
* **A memoised function obeys the generator contract** (the user's ruling
  of 2026-10-04, "Contract + cold-rebuild control"): it reads data only by
  opening real files, depends on a file's presence and bytes only, reads
  no environment variable, starts no process, branches on no other memo's
  ``lookup``, and takes arguments with a constructor form. A recorder
  written in Python cannot be airtight against arbitrary Python; the
  contract states what it does not try to catch, and phase P5's cold
  rebuild is its witness.
* **What no manifest can see is stated, not gated**: what the contract
  excludes, native state, and a file a C library opens itself. Every rule
  here errs toward a spurious miss, never toward a stale hit, except at
  those named limits.


What the memo is for
====================

Phase P3 of #405 has one job, set by the phase line and by ruling 3 of
2026-10-03 ("certification is decoupled from caching",
:ref:`verification-reference-architecture`): cache every reference
reading, certified or not, under a key that misses exactly when something
that ran to produce it has changed. The user named the reason: the
development cadence. A reference reading is a pure function of its
question and of the code that answers it, so recomputing it on every test
run buys nothing as long as neither has changed, and it costs most of the
wall time of the verification suite.

"Exactly" has two halves, and the design is held to both:

* **No stale hit.** If anything that ran to produce the answer has
  changed, the entry is not served. This half is correctness: a stale hit
  would make a test pass against the answer of code that no longer exists.
* **Few spurious misses.** If nothing that ran has changed, the entry is
  served. This half is the cadence: an entry invalidated by every edit to
  the package is a cache that never hits.

Where the two conflict, the design takes the miss. Every declared
spurious miss below (an unread module constant, a distribution upgrade
that changed nothing the answer used, an operating-system update) is
priced and accepted; the only stale hits are the named limits of "What no
manifest can see".


.. _verification-reference-cache-identity:

Why a recorded execution: the identity of a generator
=====================================================

The question of discussion 2 of the plan (opened 2026-09-24): what must
change for a cached reference to be stale? Four candidates were
considered.

.. list-table:: The four candidate identities of a generator
   :header-rows: 1
   :widths: 22 40 38

   * - Candidate
     - What it is
     - Verdict, with the structural reason
   * - A hand-declared version
     - A version string the generator's author bumps when the answer
       changes.
     - Refuted. It is prose standing in for enforcement
       (``instrument-doctrine`` X3): nothing fails when an author forgets
       to bump it, so it goes stale in silence.
   * - A module-level source-closure hash
     - The hash of every first-party file the generator's module imports,
       transitively.
     - Refuted as too coarse. `[M]` 2026-09-24, a static AST closure over
       the 353 tracked ``orpheus/**/*.py``: the cylinder multi-region
       Green's function module's closure is 124 files (59 430 lines); the
       Peierls geometry module's 114, the S\ :sub:`N` face-transmission
       module's 116, the Sood registry's 113. Every closure is about a
       third of the package and nearly the same third, because the package
       ``__init__`` files import eagerly (``orpheus/derivations/__init__.py``
       imports ``reference_values``, whose closure alone is 112 files). An
       edit to any of about 88 shared files in ``numerics``, ``geometry``
       and ``data`` would invalidate every reference: narrower than the
       git commit, and in the same failure class.
   * - A static function-level call graph
     - The functions reachable from the generator by the call graph.
     - Refuted as unsound. The project's own graph tool measured 12 to
       15 % recall of executed call edges (properties, dunders, callbacks
       and polymorphic dispatch mint no static edge), so a missed
       dependency would serve a stale reference in silence, which is the
       one failure the cache must never have.
   * - The recorded execution
     - The functions that actually ran, with the files read and the
       distributions used, recorded while the generator runs.
     - **Ruled**, 2026-09-24 (the user: "Option 4 is excellent"). It is
       as fine as the run itself: a function that did not run cannot have
       changed the answer, and a function that ran is pinned whatever
       route reached it.

**Recording costs nothing measurable.** `[M]` 2026-10-04 (``main``
``588646f0``, one run each): one eigenvalue reading of the A|B|A sphere at
``n_r=8, n_mu=8, n_traj_quad=16`` took 3.97 s untraced and 3.52 s traced by
``sys.monitoring`` with a ``PY_START`` callback that disables itself after
the first report of each code object; the difference is noise. Hashing
every recorded code object took 0.01 s. The trace held 1 053 code objects:
192 first-party (``derivations`` 92, ``numerics`` 45, ``data`` 26,
``geometry`` 16, ``specification`` 8, ``reference`` 5), 672 from the
standard library or third-party distributions, and 189 generated
(``<string>``, the dataclass-generated ``__init__`` and ``__eq__``, and the
frozen import machinery).

**A recording sees only code that runs**, and that is the fact that
shaped the rest of the design. Three things serve an answer from memory
without running the code that produced it:

1. **A function-cache hit.** A ``functools.cache`` or ``lru_cache`` hit
   inside a generation returns without running the cached function, so the
   record omits it and a later edit to it would serve a stale entry.
2. **A** ``cached_property`` **computed earlier**, on an object that
   outlives the generation that computed it (an object returned by a
   function cache, or held at module level): the same omission.
3. **A monkeypatch.** The record holds the spy's code, and the spy's entry
   validates against the checkout for as long as the test file is
   unchanged; a later unpatched call with the same key would then be served
   the spy's answer.

`[M]` 2026-10-04, by AST over ``orpheus/``: 26 decorated function caches
(``functools.cache`` 16, ``lru_cache`` 10; hand-rolled dictionary memos not
counted) and 101 ``cached_property`` sites, and 10
``monkeypatch.setattr`` in the reference gates, among them a spy that
replaces both multi-region solvers.


.. _verification-reference-cache-fresh-process:

Why a fresh process for every miss
==================================

Two architectures answered the three hazards.

**A, an in-process memo traced from session start.** One traced-memo
primitive would replace every in-process cache; a hit would replay its recorded trace into every active trace; a
gate would ban ``functools`` caches (and possibly ``cached_property``) in
``orpheus/``, and 26 to 127 sites would migrate; a separate guard would
refuse to write an entry while any traced function is not the one bound in
its module.

**B, one fresh process per miss.** A reading is looked up in the process
that asks. On a miss it is generated in a new interpreter that records
from its first line, writes the entry and exits.

**The user ruled B** (2026-10-04, question 1): *"Architecture B: one fresh
process per miss."* The reasons:

* **Hazards 1 and 2 cannot arise.** A process born for one call has
  nothing in memory to serve from, so every function the answer depends on
  runs inside the recording. This holds by construction, not by a gate
  over 127 sites, and it holds however production code chooses to memoise
  later. There is accordingly no gate on ``functools`` caches or
  ``cached_property`` in ``orpheus/``: the fact such a gate would rest on
  (an in-process cache hit hides its dependency from a trace) is true, and
  a fresh process makes it harmless.
* **Hazard 3 is safe by default.** A monkeypatch in the test process does
  not reach the generating process, which imports the code from disk, so a
  patched function can never write an entry. The test that patches must
  instead read without the cache (the section "The bypass" below), which
  is explicit.
* **The cost is paid only on a miss.** `[M]` 2026-10-04, ``-O``, three
  runs: importing the trajectory-resolvent reference, the verification
  module and the specification module in a new interpreter took 0.92,
  0.65 and 0.62 s, against generation times of seconds to minutes. Once
  the cache is warm, misses are rare, so architecture A's only advantage,
  no process boundary, rarely matters.

**What crosses into the process is content, never a built object.** The
census of 2026-10-04 (finding F1) measured the trap: unpickling a built
``TrajectoryResolventDerivation`` or ``Billiard`` restores the object
without running its construction, so in 2 of 2 traced trajectory runs the
construction functions (``Billiard.__post_init__``, ``_route``,
``_layered_xs_payload``, ``_read_isotropically``, ``reference_body`` and
the derivation's ``__post_init__``) were absent from the record; in the
positive control, which constructs the derivation inside the process, 7 of
7 were present. So every dataclass crosses through its constructor
(``_ConstructorPickler`` rebuilds it from
:func:`~orpheus.numerics.content.constructor_arguments`, and a
:class:`~orpheus.numerics.content.ContentIdentity` value already pickles
that way), and its construction runs, and is recorded, in the generating
process. The pickle is safe in this direction: the parent wrote it and the
child is our own process. The pickle ban applies to what is read back from
the cache directory, which continuous integration may restore from another
machine.

**A miss inside a generation recurses.** A memoised call made by a
generating process is looked up like any other: a hit is served and
pinned as a child of the entry being generated; a miss starts its own
generating process. So every entry carries exactly its own record, and
the solve that several readings share is one entry, pinned by each reading
that read it.


.. _verification-reference-cache-key:

The key: exact, not content equality
====================================

The key is the 256-bit BLAKE2b digest of

.. code-block:: text

   encode_exact((function id, signature-bound arguments, platform tag))

where the function id is ``module:qualname``.

**The arguments are bound to the signature**, defaults applied
(``inspect.signature(f).bind(...).apply_defaults()``). Without the binding
one call has two keys: ``Billiard`` omits the solver settings it leaves at
``None``, while a direct caller may pass the same default explicitly
(census, 2026-10-04).

**A parameter may declare a canonical form**, applied before the key is
taken, and the generating process receives the canonical argument, not
only the key: ``traced_memo(canonical={"radii": as_float_array})``. The two
multi-region solvers declare their six array parameters (``radii``,
``sigma_t``, ``sigma_s``, ``nu_sigma_f``, ``chi``, ``initial_psi``) in
``SOLVER_ARRAY_ARGUMENTS`` as ``np.asarray(·, dtype=float)``, which is the
form each solver's body first reads them in. A caller passing a list of
radii and a caller passing the array therefore share one entry, and the
answer is unchanged because the solver would have converted it anyway. A
canonical key naming no parameter of the function is refused at
decoration. The alternative, normalising every array-like globally, was
rejected because it is unsound for a parameter where a list and a tuple
mean different things.

**Exact against content equality** (the user's ruling of 2026-10-04,
Q2). The content encoder :func:`~orpheus.numerics.content.encode` follows
``==`` by the user's ruling of 2026-10-02
(:ref:`structured-geometry-content-identity`): two values that are the
same physics are one value. A function can still tell such values apart,
so a key that followed ``==`` could serve one call another call's answer.
:func:`~orpheus.numerics.content.encode_exact` is the same walk of the
value (an object is still its schema and its parts, a mapping is still
ordered by its encoded keys) with an exact leaf policy: every scalar with
its type and its bits, every array with its dtype and its bytes.

.. list-table:: Content equality and the memo key on nine pairs
   (`[M]` 2026-10-04, ``encode`` and ``encode_exact`` at ``994ba740``)
   :header-rows: 1
   :widths: 30 20 20 30

   * - Pair
     - ``encode`` equal
     - ``encode_exact`` equal
     - What a function can see
   * - ``8`` and ``8.0``
     - yes
     - no
     - ``type(x)``; ``x // 3`` is ``2`` or ``2.0``, and the payload
       returns the type computed
   * - ``-0.0`` and ``0.0``
     - yes
     - no
     - ``math.copysign(1, x)``
   * - ``np.arange(3)`` and its float copy
     - yes
     - no
     - the dtype, and so integer division and overflow
   * - ``True`` and ``1``
     - yes
     - no
     - ``type(x)``, ``repr``
   * - ``S(1.0, tag="a")`` and ``S(1.0, tag="b")``, ``tag`` a
       ``compare=False`` field
     - yes
     - no
     - ``s.tag``
   * - a ``NamedTuple`` ``NT(1, 2)`` and ``(1, 2)``
     - yes
     - no
     - ``t.a``, ``type(t)``
   * - ``FrozenMapping`` of ``a, b`` and of ``b, a``
     - yes
     - no
     - iteration order
   * - a CSR matrix and its twin storing an explicit zero
     - yes
     - no
     - ``m.data``, ``m.nnz``
   * - a CSR matrix and its CSC copy
     - yes
     - no
     - ``m.format``, ``m.indices``

**Exact all the way up.** The exactness must hold at every depth and in
every structure, not only at the leaves: a call on ``Spec(8.0)`` must not
be served the entry of ``Spec(8)`` when ``Spec`` is a frozen dataclass,
and ``sign(-0.0)`` must not be served the entry of ``sign(0.0)``. So the
exact policy rides the content encoder's own walk, one walk with two
policies (``_Exactly`` beside ``_ByEquality`` in
:mod:`orpheus.numerics.content`), rather than a second, shallower
description of the arguments (review finding qa 4). The exact policy
decides more than the leaves:

* **an object is its constructor form**,
  :func:`~orpheus.numerics.content.constructor_arguments` (or the
  arguments of its ``__reduce__``), which is exactly what the generating
  process receives, since every argument crosses through its
  constructor. A field declared ``compare=False`` that the constructor
  takes is therefore in the key, although content equality ignores it;
* **a container carries its concrete type**: a ``NamedTuple`` is not the
  plain tuple of its values;
* **a mapping keeps its order**:
  :class:`~orpheus.numerics.content.FrozenMapping` documents its iteration
  order as behaviour, so a function can read it;
* **a sparse matrix is its format and its stored arrays**, explicit zeros
  included;
* **a subclass of** ``str``, ``bytes`` **or** ``ndarray`` **is its type and
  its own instance state beside its bytes**, the state including
  ``__slots__``. ``EmissionSpectrum``, the stateless ``ndarray`` subclass
  ``Mixture.chi`` holds, is its type and its bytes; a masked array keeps a
  type in its state, which has no content identity, and is refused with
  its path;
* **a subclass of** ``int`` **or** ``float`` **is its type and its state
  beside its value**, the value taken through ``int(value)`` or the float's
  bits, never through ``__str__``, which a subclass can override;
* **a set keeps its iteration order**: two equal sets can iterate
  differently, and a function that iterates one can tell them apart;
* **a sparse matrix is keyed by its class as well as its format** (a
  ``csr_matrix`` is not a ``csr_array``), and one that stores no arrays (a
  ``dok_matrix``) is refused;
* **a mapping other than** ``dict`` **and**
  :class:`~orpheus.numerics.content.FrozenMapping` **is refused**: a
  ``defaultdict`` holds a factory beside its items;
* **a dataclass with an** ``InitVar`` **is refused**: its constructor
  takes a value it does not keep, so it has no constructor form.
  :func:`~orpheus.numerics.content.constructor_arguments` refuses it, which
  also refuses it at the process boundary and in the payload.

The content policy is unchanged by all of this: every content digest is
bit-identical (`[M]` 2026-10-04, the S5.7 RECORD fingerprint and the
content-identity gates green: 244 at ``994ba740``, 224 named in
``bacf0787``'s message). The price is a spurious miss for a caller
passing ``8`` where another passed ``8.0``, or one dictionary order where
another passed the other, which is the safe direction.

**The platform tag** (Q3) is
:func:`~orpheus.numerics.traced_memo.platform_tag`: the interpreter's
cache tag, the operating system, the machine, the ``-O`` level and the
name of numpy's BLAS. The ``-O`` level is in it because a bare ``assert``
in a generator is stripped under ``-O``, so the two levels can answer
differently; the BLAS is in it because Apple Accelerate and OpenBLAS
differ in the last bits (#504). The BLAS thread count is not in it: whether
it moves the last bits here is unmeasured, and the gate that would show it
is M4.5 (a served reading equal bit for bit to a fresh in-process
reading, on one machine).

**A memoised method** is keyed on its instance. ``TracedMemo.__get__``
binds a memo like a function, so ``derivation.evaluate(observable)`` is
the memo called with ``(derivation, observable)`` and the derivation's
content is part of the key; the instance must therefore have content
identity. An argument with no content identity is refused with its path
(M3.6), never keyed by its address.


.. _verification-reference-cache-manifest:

The manifest: eleven kinds of pin
=================================

A :class:`~orpheus.numerics.traced_memo.Manifest` is a frozen set of
pins. Each pin is a small value that can say, by itself, whether the
current checkout still matches what it recorded (its ``stale`` method
returns a reason or nothing). The JSON form tags every row with its pin's
kind and is parsed once, in ``Manifest.from_json``.

.. list-table:: The pin kinds
   :header-rows: 1
   :widths: 17 33 25 25

   * - Pin
     - What it pins, and by what
     - What it guards against
     - Where it came from
   * - :class:`~orpheus.numerics.traced_memo.DefPin`
     - A first-party def that ran: its file, its qualified name, its
       rank among the defs of that name, and the digest of its normalised
       source.
     - An edit to the body of a function that ran.
     - The design; the rank from qa finding 1.
   * - :class:`~orpheus.numerics.traced_memo.ModulePin`
     - A first-party module whose code ran, by its skeleton.
     - An edit to a constant, an import, a decorator, a signature default
       or a class attribute that a traced function reads by name.
     - The design.
   * - :class:`~orpheus.numerics.traced_memo.DistributionPin`
     - A third-party distribution whose code ran or whose file was read,
       by its installed version.
     - A dependency upgrade (numpy, scipy, mpmath, SymPy).
     - The design.
   * - :class:`~orpheus.numerics.traced_memo.InterpreterPin`
     - The interpreter's version and bytecode tag.
     - The standard library and generated code changing under a new
       Python.
     - The design.
   * - :class:`~orpheus.numerics.traced_memo.DataPin`
     - A file the run opened for reading, by absolute path and bytes;
       ``"absent"`` when it was missing; a ``.py`` file read as data
       included. Also an admitted program the run started (``uname``), by
       its executable's bytes.
     - An edited data file; a file appearing where the run found none; a
       replaced program.
     - Q4 (M3.17); absence from qa finding 5; ``.py`` as data from the
       second review's finding r5; programs from the fix rounds.
   * - :class:`~orpheus.numerics.traced_memo.PresencePin`
     - A path the run probed without opening it (``os.stat``,
       ``os.lstat``, and through them ``os.path.exists``, ``isfile``,
       ``isdir`` and ``pathlib``'s ``exists``), by what was there:
       ``"file"``, ``"directory"``, ``"other"`` or ``"absent"``.
     - An answer that branched on whether a file exists.
     - qa finding 5, the half the first fix left open.
   * - :class:`~orpheus.numerics.traced_memo.ListingPin`
     - A directory the run listed, outside the import system, by the
       sorted names in it (or ``"absent"``).
     - A file added to or removed from a directory a generator scanned.
     - qa finding 5.
   * - :class:`~orpheus.numerics.traced_memo.EnvironmentPin`
     - The declared variables the generating process received, ``PATH``
       excluded, by their digest.
     - A thread count or a locale that changed what ran.
     - The second review (an undeclared variable read by a generator).
   * - :class:`~orpheus.numerics.traced_memo.WorkingDirectoryPin`
     - The directory the run started in, recorded only when the run read a
       relative path.
     - The same relative path naming another file from another directory.
     - qa finding 5; the starting directory from finding r7.
   * - :class:`~orpheus.numerics.traced_memo.ChildPin`
     - A memo entry the run read, by its function id, its key and the
       digest of the payload it served.
     - A child entry that went stale, or was regenerated with another
       answer.
     - The design (recursive validation).
   * - :class:`~orpheus.numerics.traced_memo.MemoPin`
     - The memo's own two source files and ``content.py``, whose
       constructor form the write path runs.
     - An entry written by a memo whose writing code has since been
       fixed.
     - qa finding 8, both halves.

**The normalised source of a def.** A def is hashed as the dump of its
syntax tree with line and column attributes excluded and every docstring
removed, so a comment, a layout change or a docstring edit never moves the
digest (the user's ruling of 2026-10-04, question 4: a docstring edit
cannot change an answer). Hashing source rather than bytecode lets
validation parse the current files without importing them.

**One walk finds every outermost def.** The pinned unit is the outermost
def: a nested function, a lambda or a comprehension inside a def is pinned
by its enclosing def. The walk descends every statement except a def's
body, so a def under any module-level block (``if``, ``for``, ``while``,
``match``, ``try``, ``with``) is found, and a class adds its name to the
qualified names of the defs in its body. The walk names no block type: a
list of the blocks to descend is a list that can miss one, and a def it
misses is never pinned, so an edit to it serves a stale hit (review
finding qa 6).

**A def is pinned by its rank among defs of one name.** Several defs can
share a qualified name: a property's getter and setter, the two arms of a
module-level ``if``/``else``, ``functools.singledispatch`` registrations
named ``_``. A lookup by name alone finds the LAST def of that name, the
one Python binds, which need not be the one that ran; an edit to the one
that ran would then serve a stale hit (review finding qa 1). So the pin
carries the def's rank among the defs of its name (0 for the first), and
the qualified name with the rank is the pin.

**The skeleton of a module** is its syntax tree with every function body
replaced by ``pass`` and every docstring removed: its imports, constants,
decorators, signatures with their annotations and defaults, class
attributes and dataclass fields. A function reads ``_EXACT_DIGITS = 60``
by name, so its own source does not change when the constant does, and the
module's skeleton does. The skeleton is pinned per module, not per name
read, so an edit to a constant no traced function reads is a **declared
spurious miss** (the row ``constant-unread`` of M2.1): sound, and priced.
The skeleton of every first-party module whose body ran enters the
manifest; the census counted 103 to 121 such modules per generation,
because the generating process imports from its first line.

**Code that is not a def body.** Python 3.14 generates code objects from
a signature or a class body rather than from a def body: lazy annotations
(PEP 649, ``__annotate__``, run whenever anything reads a function's
annotations) and the scopes of type parameters (PEP 695,
``<generic parameters of …>``). Their source is the signature, which the
skeleton pins, so they never pin the body of the def they belong to (M1.3c):

* an ``__annotate__`` is reported at its def's line, so mapping every code
  object to the def spanning its line would pin a def whose body never
  ran whenever anything read its annotations (the payload decoder's
  ``typing.get_type_hints`` does);
* a type-parameter scope such as ``<generic parameters of find_factor>``
  (``orpheus/numerics/space.py``) has no def of its own name at its line,
  so a pin that demands one would refuse every generation that imports
  ``space.py``.

Code compiled under a file's name at a line where the file holds no such
code is refused as :class:`~orpheus.numerics.traced_memo.Unpinnable`
(M1.6): it ran, and no pin can describe it.

**Where a file belongs**, decided once by
:func:`~orpheus.numerics.traced_memo.origin`: generated code
(``<string>``, ``<frozen …>``) is dropped, because the interpreter pin and
the skeleton that generated it cover it; a standard-library file is
covered by the interpreter pin; a file in ``site-packages`` is pinned by
its distribution's version, found through
``importlib.metadata.packages_distributions`` and, for a file no top-level
name maps, through the distributions' ``RECORD`` files; everything else is
source, pinned by content. A file of an **editable** distribution is pinned
as source, because an editable install's version never moves: every
generating process runs the editable ``orpheus`` finder,
``__editable___orpheus_0_1_0_finder.py``, which no top-level name maps
and which only the ``RECORD`` search places. Searching the raw ``RECORD``
texts for one file costs 0.23 s; building the whole file-to-distribution
index through ``importlib.metadata`` costs 9.4 to 10.9 s per process
(`[M]` 2026-10-04, 24 484 files). A ``site-packages``
file that belongs to no distribution is refused as unpinnable.

**A source file is pinned relative to the** ``sys.path`` **entry that
holds it**, so validation finds the file the same import would load now;
a source file on no ``sys.path`` entry (one imported through the editable
install's finder) is pinned by its absolute path, so the memo works from
a process that imports ``orpheus`` through the finder (review finding
qa 12).

**Data files** (Q4, gate M3.17). The trace records code, not
the files a generation reads. An audit hook on ``open`` in the generating
process records every path opened for reading, made absolute at the moment
it is opened (a later change of directory cannot move it), a bytes path
decoded, and a path that does not exist included: the run then depended on
its absence. A ``.py`` file read as data (``inspect.getsource`` of code
that never ran, a parsed configuration) is a data file like any other;
only the import system's own reads are exempt (the next paragraph).

**Presence probes.** While a recording is active, ``os.stat`` and
``os.lstat`` are replaced by probes that record the path, so
``os.path.exists`` and its relatives are dependencies: each probed path
becomes a :class:`~orpheus.numerics.traced_memo.PresencePin` holding what
was there. The probes are removed when the recording stops. In P3, 0 of the 3 clients read a nuclear-data file (`[M]` the
AST census: 4 files under ``orpheus/`` read data, all in ``micro_xs`` and
the Sood cache), so the hook is for P4's generators.

**Directory listings** (``os.listdir``, ``os.scandir``) are recorded.

**The import system is exempt, decided by the immediate caller.** The
import system's own reads, probes and listings are how it finds a module,
and the module it finds is pinned where its code runs; recording them
would let any new file in any ``sys.path`` directory invalidate every entry
(`[M]` 2026-10-04: the M4.4 real-tree witnesses fail that way without the
exemption). The import system is the frozen bootstrap and the
``importlib`` package (``importlib.metadata`` lists every ``sys.path``
directory to find the installed distributions, which are pinned by
version). The exemption is decided by the frame that made the call, never
by an import being somewhere on the stack: a module body that lists a
directory while it is being imported is the run's own dependency (finding
r4). The check compares raw file names against prefixes computed once,
because it runs inside the ``os.stat`` probe and must not call
``os.stat`` itself.

**The memo's own store reads are not dependencies.** A memoised call made
inside a generation looks its entry up, and that lookup validates the
entry's manifest by reading each of its source files. The lookup runs
quiet (``_traced_memo_boot.unrecorded``, a thread-local flag the audit
hook and the ``os.stat`` probe both honour), so a parent holds its child
by its :class:`~orpheus.numerics.traced_memo.ChildPin` only, and an edit
to one of the child's source files reaches the parent through the child's
own validation, never as a data pin of the parent.

**Programs a run starts.** ``platform.processor()`` runs ``uname`` during
the clients' imports (`[M]` 2026-10-04), so a started program is a real
dependency. A generation may start only a declared machine query, a
program whose answer the machine fixes: ``_ADMITTED_PROGRAMS`` holds
``uname`` by its system path (``/usr/bin/uname`` or ``/bin/uname``,
resolved), never by its name, so a script named ``uname`` earlier on
``PATH`` is refused (third review); it is pinned by its executable's bytes
(an operating-system update can replace it). Every other start is refused as
:class:`~orpheus.numerics.traced_memo.Unpinnable`, because a started
program's own reads and children are unrecorded: a pinned ``cat`` of an
unpinned file served ``1.0`` where ``50.0`` was computed, and
``shell=True`` and ``/usr/bin/env python3`` slipped past a refusal keyed on
the program being Python (findings r2, r3). A spawn-context process pool,
which raises the audit event ``_posixsubprocess.fork_exec`` rather than
``subprocess.Popen`` and is named by its executable candidates, a ``fork``
and an ``os.system`` line are refused the
same way. The table is a declared scope boundary (``SCOPE-BOUNDARY`` in
the code): the machinery that would admit more is a recorder for a
started program's own reads and children. The memo's own start of a
generating process is not a dependency of the run that asked, and is
excluded by name (``_traced_memo_boot.spawning``).

**The memo pins itself.** The memo writes an entry after the recording
stops, so its own writing code is never in the record. Without a pin on
its own source, an entry written by a memo with a codec defect would stay
a ``Hit`` after the defect was fixed (review finding qa 8). The pin covers
the memo's two source files and ``content.py``, whose
:func:`~orpheus.numerics.content.constructor_arguments` the write path
runs; an edit to any of them makes every entry stale.

**What a real manifest holds** (`[M]` 2026-10-04 at ``6bfde984``, this
chapter's probe, the A|B|A sphere at ``n_r=8, n_mu=8, n_traj_quad=16``,
eigenvalue reading, fresh cache root; every entry also holds one
interpreter pin, one memo pin and one environment pin):

.. list-table::
   :header-rows: 1
   :widths: 36 9 9 9 9 9 9 10

   * - Entry
     - Def
     - Module
     - Distr.
     - Data
     - Presence
     - Child
     - bytes
   * - the reading (``TrajectoryResolventDerivation.evaluate``)
     - 208
     - 124
     - 6
     - 4
     - 55
     - 1
     - 51 147
   * - the solve (``solve_greens_function_sphere_mr``)
     - 62
     - 105
     - 6
     - 4
     - 55
     - 0
     - 31 335

The four data pins of both entries are
``/System/Library/CoreServices/SystemVersion.plist`` (read by a platform
query during the import closure), ``/dev/null``, the interpreter's
``python314.zip`` (an entry of ``sys.path`` that does not exist, so
pinned as ``"absent"``) and ``/usr/bin/uname`` (the program ``platform.processor``
starts). An operating-system update therefore misses every entry; that is
sound, since such an update can change ``libm`` and the BLAS (the
specification's finding O). The six distributions are what the import
closure ran, not only what the answer used (``charset-normalizer``,
``h5py``, ``mpmath``, ``numpy``, ``scipy``, ``setuptools`` on the
prototype's manifest): each is a spurious miss on an upgrade, never a
stale hit.


.. _verification-reference-cache-validation:

Validation runs nothing
=======================

:func:`~orpheus.numerics.traced_memo.validate` asks every pin whether the
checkout still matches it and returns every reason it does not, sorted. It
imports nothing and runs nothing (M2.4 counts it): a def pin parses the
current file and digests the def; a module pin digests the skeleton; a
distribution pin reads the installed version's metadata; a child pin looks
the child up in the same store, recursively. The parse of a file is cached
by the file's BYTES, so an edited file is a new key and never a stale read
(M2.5), and in one process the first validation costs about 0.3 s and every
later one about 10 ms (`[M]` 2026-10-04 on the prototype: hashing the 121
first-party files of the sphere workload 8.5 ms per pass, parsing them
once 268.5 ms).

A lookup's verdict is the closed sum ``Hit | Absent | Stale | Corrupt``:

* :class:`~orpheus.numerics.traced_memo.Hit`: the entry validates. It
  carries the payload digest, the payload tree and the arrays exactly as
  they were verified, and those are what is decoded and served: a second
  read of the file after verifying it could meet an interleaved write and
  serve a value no generation produced (review finding qa 7).
* :class:`~orpheus.numerics.traced_memo.Absent`: no entry.
* :class:`~orpheus.numerics.traced_memo.Stale`: the manifest no longer
  describes the checkout, with every reason (``function <path>:<qualname>
  changed``, ``module <path> skeleton changed``, ``child …: its payload …
  became …``).
* :class:`~orpheus.numerics.traced_memo.Corrupt`: the entry cannot be read,
  or its payload does not match its digest.

Only a ``Hit`` is served; every other verdict regenerates and overwrites.
A freshly generated entry that does not validate is an error, never served.
The verdict is read through ``memo.lookup(*args, **kwargs)`` without
generating, which is how every gate observes a hit or a miss: from the
verdict and its reasons, never inferred from timing.

**Recursive validation and early cutoff.** A child pin compares the
child's current payload digest with the one the parent read. A child that
went stale makes its parent stale through the pin (M3.12b); a child
regenerated with the SAME answer leaves the parent valid, because the
digest is over the payload's content (its tree and each array's name,
dtype, shape and bytes), not over the zip container that stores it; a
child regenerated with ANOTHER answer makes the parent stale (M3.12c). In
the vocabulary of build systems the memo is a build with dynamic
dependencies (recorded while the task runs) and early cutoff (a rebuilt
dependency with an unchanged result does not rebuild its dependents).


.. _verification-reference-cache-payload:

The payload: exact, typed, never a pickle
=========================================

The answer is written by :func:`~orpheus.numerics.traced_memo.encode_payload`
and read by :func:`~orpheus.numerics.traced_memo.decode_payload`. It admits
exactly these types, and nothing that merely subclasses one of them:

* ``None``, ``bool``, ``str``, ``int`` (as its decimal text) and
  ``float`` (as ``float.hex``, so every bit returns);
* a numpy scalar of boolean, integer, real or complex kind, with its
  dtype (``np.float64`` subclasses ``float``, so a codec that tests
  ``isinstance(v, float)`` first would hand ``k_eff`` back as another
  type; M2.7);
* a numpy array of those kinds, stored by name in the entry's ``.npz``
  with its dtype, shape and bytes;
* a tuple of admitted values;
* a dataclass the function's return annotation names, transitively through
  its fields, stored by its constructor arguments and rebuilt through its
  constructor, so its laws re-run on load.

**Exact means exact type.** A codec that admits by ``isinstance`` loses
what makes a subclass the subclass: a ``NamedTuple`` comes back a
``tuple``, an ``IntEnum`` an ``int``, and a masked array loses its mask
(review finding qa 3). A subclass of an admitted type is therefore refused
with :class:`~orpheus.numerics.traced_memo.Unencodable`, never written as
its base.

**Only declared dataclasses.** A payload may name only a type the
function's return annotation names, so loading a payload never imports a
type its function does not declare. An undeclared dataclass is refused at
the WRITE, in the generating process, and the error crosses to the caller;
refused only at the load, the entry would be generated, written and refused
on every call (review finding qa 11b).

**Never a pickle.** Loading a pickle executes code, and the cache directory
may be restored by continuous integration from another run (the plan's
discussion 2, second exchange). The entry is read with
``allow_pickle=False``.

**One file, replaced atomically.** An entry is one ``.npz`` file at
``<root>/<function id>/<key>.npz``; its JSON (the function id, the key,
the manifest, the payload tree and its digest) is the member
``__entry__``. It is written to a temporary file beside it and moved into
place with one ``os.replace``, so a reader sees the old entry or the new
one and never a part of either. A file is used rather than a directory
because no single call replaces a non-empty directory atomically (review
finding qa 10).

**Every load is fresh and read-only.** Each served array is a new copy
with its write flag cleared. A consumer that wrote into a served array
would otherwise corrupt the next caller's hit, so the read-only flag is
also the detector of such a consumer.
`[M]` the specification's acceptance run (M4.7, 14 target files, 785
rows, cold and warm): no consumer wrote into a served array (a write
raises).


.. _verification-reference-cache-boundary:

The process boundary
====================

A generation is a ``_Job`` (the memo, the bound arguments, the key, the
store's root, and the chain of calls the generating processes are already
answering), pickled through the constructor pickler and sent on the
standard input of

.. code-block:: text

   <the parent's interpreter> [-O…] -P orpheus/numerics/_traced_memo_boot.py

with the parent's ``sys.path`` and working directory, and a declared
environment only (below).

**The boot script imports only the standard library**, so it can start
recording before any ``orpheus`` code runs: the first statement of its
``main`` starts a :class:`~orpheus.numerics._traced_memo_boot.Recording`,
and only then does it load the job and the memo. It is run by path with
``-P``, because the script's own directory would otherwise be first on
``sys.path``, and ``orpheus/numerics/operator.py`` would shadow the
standard library's ``operator``.

**One recorder.** :class:`~orpheus.numerics._traced_memo_boot.Recording`
serves both the generating process and
:func:`~orpheus.numerics.traced_memo.trace_call`, the in-process
instrument the manifest gates use, so the gates exercise the production
recorder (elegance finding S2). A recording holds, across every thread of
the process, every code object that started, every path opened for
reading or probed, every directory listed, every program started, the
directory the run started in, and every memo entry read. The state is process-wide, not held in context variables: a
new thread does not inherit context variables, so a memo called from a
worker thread would be neither pinned as a child nor written to the right
store (review finding qa 2). Audit hooks
cannot be removed, so one hook is installed per process and serves the
recordings active at each event; the ``os.stat`` probes are installed
while a recording is active. A process pool started by a generator is a
started program outside the admitted table and is refused. The builtin
``exec()`` raises the audit event ``exec`` (a dataclass's generated
methods are built that way), and it is not a process start, so it is not
a spawn event.

**The boot module is registered under its own name.** In the generating
process the boot file runs as ``__main__``, and the memo imports it by its
module name. Without an alias, that import would load a second copy of the
boot module whose list of active recordings is empty, so every memo entry
the run read, and the memo's own start of a nested generating process,
would be offered to nobody. The boot's ``main`` therefore registers itself
in ``sys.modules`` under its module name before anything else.

**The answer travels on a reserved descriptor.** The boot script
duplicates standard output for the answer and points standard output at
standard error before the job runs, so a ``print`` in a generator cannot
corrupt the answer (``orpheus/data/macro_xs/mixture.py`` prints
progress with ``flush=True`` twice; review finding qa 9).

**An exception crosses back.** A generator that raises writes no entry and
the caller re-raises the same exception, with a note naming the
generating process (M3.11). An exception that pickles in the child but
cannot be rebuilt in the parent crosses as a ``RuntimeError`` carrying the
child's formatted traceback, never as an anonymous failure (review
finding qa 11a).

**A memo calling itself with its own key is refused** with
``RecursionError``: each job carries the calls its ancestors are
answering, and a generating process that meets its own ``(function id,
key)`` among them raises. Without it the cycle starts interpreters without
end (review finding qa 13; `[M]` 2026-10-04, the fix battery's arm that
removes the guard left 192 recursing processes). A chain through DISTINCT
keys (``f(n)`` calling ``f(n + 1)``) repeats no key, so the depth of
nested generations is also bounded, at 16
(``_MAX_GENERATION_DEPTH``); the 17th raises ``RecursionError``.

**The environment is declared** (Q7, then the second review). The generating
process receives only ``_CHILD_ENVIRONMENT`` (``PATH``, ``HOME``,
``TMPDIR``, the locale variables and the BLAS and OpenMP thread counts,
each when the parent has it) and ``PYTHONHASHSEED=0``, fixed so that
iteration over a set of strings has one order in every generation. The
pin, an :class:`~orpheus.numerics.traced_memo.EnvironmentPin`, holds the
declared variables as the run RECEIVED them, ``PATH`` excluded: ``PATH``
is passed (``uname`` is found through it) but the one program a
generation may start is admitted by its system path and pinned by its
bytes, and a pinned ``PATH`` made every entry stale when a virtual
environment was activated; and a pin taken after the run would never
validate for a generator that sets a declared variable (third review,
``test_e_…``). Two reasons: no
variable may select an answer that no key or pin holds (a generator
reading an undeclared variable ``QA_SCALE`` served ``2.0`` where ``3.0``
was computed; ``test_m3_16_…``), and the withdrawal switch
``ORPHEUS_RUN_WITHDRAWN`` must not let a withdrawn generator write an
entry that a later caller, without the switch, would be served (the
specification's finding W). No ``ORPHEUS_*`` variable is declared, so a
withdrawn generator always refuses inside a generating process, and the
rule P4 inherits is that a withdrawn generator is never memoised.


.. _verification-reference-cache-bypass:

The bypass
==========

``with bypass():`` (:func:`~orpheus.numerics.traced_memo.bypass`, Q9)
makes every memoised call inside the block run in the calling process,
reading and writing nothing. It is the spelling for a test that
monkeypatches anything a generation runs: a patch never reaches a fresh
process, so without the bypass the test would read the honest entry
instead of its patched answer, and a spy counting calls would count
nothing. The census found two such sites: the solve spy of
``tests/gates/derivations/test_trajectory_resolvent_reference.py``, which
patches both multi-region solvers to count lazy solves (its rows read
under an ``in_process`` fixture that enters the bypass), and a counting
spy on ``Symbolic.steps``. M3.10 gates both halves: a patch in the parent
never reaches an entry, and a call under the bypass reads and writes
nothing.

The bypass is a context manager, not an environment variable, because an
environment variable is exactly the ambient state the ambient-state step
removed (below). It and
:func:`~orpheus.numerics.traced_memo.cache_root` (read and write under
another root) are entries on one process-wide stack of targets; a stack
in context variables would not reach a worker thread (review finding qa 2,
elegance finding C3). The other side of that choice is a declared limit:
a ``bypass()`` entered in one thread holds for every thread of the
process (found by the second review and declared, not fixed; low, since
the tests that bypass are single-threaded).


.. _verification-reference-cache-certificate:

The memo and the certificate
============================

The user asked, on the specification's question Q5: *"It can execute if
it's fast, but how is the certificate handled?"* The answer, from
:mod:`orpheus.reference.solution`:

* **The memo sits below the certificate.** It caches what
  ``derivation.evaluate(observable)`` returns: an establishment
  (:class:`~orpheus.reference.certificate.Exact` or
  :class:`~orpheus.reference.certificate.DerivedBound`, with its enclosure)
  or an :class:`~orpheus.reference.reading.Uncertified` value. It never
  caches a certificate, a certificate's state or a verification verdict.
* **The certificate is built in the process that asks.** The factory
  builds the :class:`~orpheus.reference.solution.ReferenceSolution` there,
  and ``ReferenceSolution.__post_init__`` re-checks every claimed
  observable against the derivation's evaluation, cached or computed: a
  claim whose evaluation leaves it uncertified, or whose established
  enclosure shares no point with the claimed one, cannot be built. A cached
  establishment that disagrees with its claim is therefore refused at
  construction, exactly as a computed one would be.
* **A reference certificate's state is derived on every read**
  (:ref:`verification-reference-architecture`), never stored, so there is
  no certificate state for the memo to serve stale.
* **The certification run of a P4 family will be a client** (the plan's
  W5 H: a certificate is a traced memo of its certification run), so the
  expensive part of certifying is memoised under the same rules, while the
  claim is still re-checked where it is consumed. `[R]`: nothing in P3
  certifies a family by a run; the one certified family, the exact
  infinite medium, is certified by exact rationals.

A wrong reference served from the cache is exactly as wrong as computed.
The memo certifies nothing about a reading, and its gates say nothing
about whether a reading is correct (the specification's section 6).


.. _verification-reference-cache-exact-medium:

Why the exact infinite medium is not memoised
=============================================

The user's ruling 3 of 2026-10-04 made the exact infinite medium a client,
for uniformity: one primitive for every producer of a
``ReferenceSolution``. The specification measured what that buys
(`[M]` 2026-10-04 on the prototype, one run each, concurrent with other
runs, so upper bounds):

.. list-table:: The exact infinite medium (mixture A, 2 groups), one
   construction and one reading
   :header-rows: 1
   :widths: 40 30 30

   * - Route
     - Time
     - Interpreters started
   * - Computed in the asking process
     - 0.01 s
     - 0
   * - Memoised, warm
     - 0.11 to 0.12 s
     - 0
   * - Memoised, cold
     - 6.30 s
     - 3
   * - ``tests/gates/reference/`` (386 rows): today / cold / warm
     - 7.6 s / 60.5 s / 14.5 s
     - —

A warm hit cost about ten times the computation, because validating an
entry parses and hashes every first-party file the generation ran, and a
cold miss about six hundred times. Memoising a computation cheaper than
its own validation buys nothing the cache exists for, and the uniformity
is already provided by the shared ``Derivation`` protocol. **The user
ruled it out** (Q5): *"It can execute if it's fast"*. The exact medium is
computed in the process that asks and is not memoised; M4.9 is the gate
that its ``evaluate`` is a plain method, starts no interpreter and writes
no entry. The certificate is unaffected either way, by the previous
section. One consequence: the specification's Q8 (the exact derivation
should hold the mixture so that its content can key a memo) does not
arise, and ``ExactInfiniteMediumDerivation`` keeps its one field,
``medium`` (`[M]` the live class).

The rule this ruling sets for later clients: **a reading cheaper than one
interpreter start plus one validation is not memoised.**


.. _verification-reference-cache-clients:

The clients
===========

**The two multi-region solvers**,
``solve_greens_function_sphere_mr`` (in
:mod:`orpheus.derivations.continuous.trajectory_resolvent.greens_function`)
and ``solve_greens_function_cylinder_mr`` (in
``greens_function_cylinder``), are decorated
``@traced_memo(canonical=SOLVER_ARRAY_ARGUMENTS)``. The census showed why
the memo sits at the solver functions themselves: they return frozen
results of 11 and 10 plain-data fields, each field read by at least one of
the 21 consumer sites (2 in production, 19 in tests); the reference's
eigenvalue and a direct call's are bit-identical (``1.1703350703390714``
on the census's fixture, `[M]` 2026-10-04);
so one entry serves the solver's direct callers and the reference, which
reaches the solver through ``Billiard``. The entry's payload carries the
knots the reading needs to rebuild its rays. This is the user's ruling 2
(the solve is a child entry shared by every observable of that solve) and
ruling 3's legacy client (the cylinder multi-region Green's function, the
first family of P4) in one decoration.

**The trajectory-resolvent reading.**
:class:`~orpheus.derivations.continuous.trajectory_resolvent.reference.TrajectoryResolventDerivation`
has content identity (its six init fields; its ``Billiard`` and its rays
class are derived from them and excluded), and its ``evaluate`` is a
memoised METHOD. A reading is generated once, in a fresh process that
constructs the derivation from its content (finding F1), solves through
the solver's memo (a child entry) and evaluates the observable. The memo
is on ``evaluate`` itself, not on a module-level reading function beside
it, which would be a second spelling of ``evaluate`` free to drift from
it. Within one process the derivation still holds its solve
in a ``cached_property``; across processes, the reading and the solve are
memo entries (:ref:`trajectory-resolvent-reference-reading`).

**A solve that does not converge.** The solver's result is a value
whatever its convergence, so an exhausted solve is a cached value; the
reading of it is a refusal, and a refusal is an exception, which writes no
entry (M4.8).


.. _verification-reference-cache-contract:

The generator contract
======================

**The user's ruling** (2026-10-04, after qa's third review), verbatim:
*"Contract + cold-rebuild control"*. The memo records what a generation
does from inside Python, and a recorder written in Python cannot be
airtight against arbitrary Python: each of three qa rounds found new ways
a deliberately adversarial generator could hide a dependency from it, and
none of those ways touched the three real clients. Hardening without end
against generators nobody writes was refused; instead the memo states the
contract a memoised function obeys, and what the contract excludes the
recorder does not try to see. The contract, as it stands in
:mod:`orpheus.numerics.traced_memo`'s module docstring:

* reads data only by opening real files (``open``, ``Path.read_*``,
  ``np.load``), never through a loader (``pkgutil.get_data``) or a
  module's absence (``try: import accel``, ``find_spec(...) is None``);
* depends on a file's presence and its bytes, never on its other metadata
  (its size or time from ``stat``, whether it is a link);
* reads no environment variable (the M5.1 census holds this for
  ``orpheus/``), starts no process, and branches on no other memo's
  ``lookup``;
* takes arguments with a constructor form: no ``InitVar``, no mapping
  beside ``dict`` and ``FrozenMapping``, no sparse matrix without stored
  arrays (each refused when keyed).

The fourth clause is enforced: each excluded argument is refused when the
key is taken. The first three are a contract, not a gate. The contract is
a declared scope boundary (``SCOPE-BOUNDARY`` in the module docstring):
the machinery that would replace it is a recorder below Python, such as a
system-call tracer. Its witness is phase P5's cold rebuild, which
regenerates every entry without the cache and compares it with the cached
one: a generator that breaks the contract and was served a stale entry
shows up there as a mismatch. A generator that needs what the contract
excludes, or a cold rebuild that finds a stale entry, reopens the ruling.

**What the contract covers** (each a way qa's third review hid a
dependency, none in a real client):

* **File metadata.** A presence pin keeps a file's kind (file, directory,
  other, absent), not its size, its times or whether it is a link, so a
  generator that branches on ``os.stat(path).st_size`` is outside the
  contract. ``os.access`` is not wrapped at all, and the entries
  ``os.scandir`` yields are not probed (the listing itself is pinned).
* **Loader reads.** A data read through ``pkgutil.get_data`` runs behind
  import-system frames, which the recorder exempts.
* **An optional import's absence.** A generator whose answer depends on
  whether an optional module can be imported depends on something the
  import system decided, and the import system's probes are exempt.


What no manifest can see
========================

Stated where the memo is built (its module docstring), and here, beside
what the generator contract excludes (the section above):

* **Native state.** A C extension's global state (a BLAS thread pool, a C
  library's precision mode) is neither code the trace records nor an
  argument. M4.5 compares a served reading with a fresh in-process one on
  one machine, and cannot see a difference both share.
* **A file a C library opens itself** (HDF5 through ``h5py``) raises no
  Python audit event. A generator reading such a file takes the file's
  digest as an argument; P4 inherits that rule.
* **A started program's own reads and children.** Only ``uname`` may be
  started; every other start is refused, not pinned.
* **A** ``bypass()`` **is process-wide**, across threads (the bypass
  section).
* **The atomic write has no gate that can fail on it.** M3.15 (two
  processes racing on one miss leave one valid entry) cannot see a
  non-atomic write, because both racers recover through the ``Corrupt``
  path (the specification's battery arm B14, declared and confirmed
  blind); qa's threaded reproducer of finding 10 is the gate
  ``test_q10_concurrent_writes_of_one_key_are_atomic``.
* **The editable install's mapping.** If the finder maps ``orpheus`` to
  another checkout, the child imports from it only when ``sys.path`` does
  not already hold ``orpheus``; M3.8 shows the path wins, and the main tree
  is always on it under ``pytest``.
* **The correctness of the answer.** The memo makes a reading cheap; it
  certifies nothing (the section on the certificate above).


.. _verification-reference-cache-ambient:

The ambient-state step: no reference reads the environment or the global precision
==================================================================================

A traced memo keys a call on its arguments and pins the code that ran;
state that is neither escapes both. The census of 2026-10-04 found four
such reads and writes in ``orpheus/``:

* ``ORPHEUS_SLAB_VIA_E1``, read from the environment at import by
  ``peierls_nystrom/cases.py``, which selected the Peierls slab route (the
  native E\ :sub:`1` Nyström path or the unified polar path). A memoised
  slab reading would have been keyed without it, so flipping the variable
  would have served the other route's answer. **The route is now an
  argument**:
  ``build_two_surface_case(..., slab_route: SlabRoute = SlabRoute.UNIFIED)``,
  with :class:`~orpheus.derivations.continuous.peierls_nystrom.cases.SlabRoute`
  (:ref:`theory-peierls-slab-polar-g5-routing`). The two rows of
  ``test_peierls_multigroup.py`` on the route pin the argument, and a
  mutation dropping the ``UNIFIED`` arm reddens both.
* The withdrawal switch (``withdrawal.py``), which selects which tests run,
  not which answer is computed. It stays, and the generating process never
  sees it (its environment is declared, above).
* **Two writes of the global mpmath precision**, ``mp.dps = 30``, in the
  cylinder Wronskian identity (``cylinder_derivations.py``) and the slab
  :math:`T_{00} = P_{ss}` identity (``greens_function_slab.py``). A global
  write leaves the precision at 30 for whatever runs next in the process,
  so an answer computed later would depend on which test ran before it.
  Both are local ``mpmath.workdps`` contexts.

The gates are ``tests/gates/test_ambient_state.py`` (M5.1 to M5.4, 14
rows): two AST censuses over ``orpheus/``, each with positive controls for
every spelling it must find (``os.environ.get``, ``os.getenv``, ``from os
import environ``, an aliased ``getenv``, an aliased ``os``; ``mp.mp.dps
=``, ``mpmath.mp.prec +=``, ``setattr(mpmath.mp, "dps", …)``, and the
negative ``with mpmath.workdps(30):``); a behavioural check that each of
the two identities leaves a precision of 17 at 17 and still passes; and a
check that two fresh interpreters, with and without the variable, see the
same module-level values. The only modules that read the environment are
the withdrawal switch and the memo, which reads it only to hand the
declared variables on.


The numbers
===========

**Cold, warm and in-process** on the A|B|A sphere at ``n_r=8, n_mu=8,
n_traj_quad=16``, one eigenvalue reading, all three readings bit-identical:

.. list-table::
   :header-rows: 1
   :widths: 40 20 20 20

   * - Measurement
     - Cold (two interpreters: the reading, the solve)
     - Warm
     - In process (bypass)
   * - Commit ``330ba98a``'s smoke test, one run
     - 7.2 s
     - 0.02 s
     - 3.3 s
   * - This chapter's probe, minimum of three runs at a load average of
       about 6
     - 10.73 s
     - 0.022 s
     - 4.32 s

The cold read costs 2.2 to 2.5 times the in-process computation (two
interpreter starts, two import closures under recording, two manifests);
the warm read costs about half a percent of it.

**The suite, before and after**, from the specification (`[M]`
2026-10-04, ``0e8740ef`` and the prototype, one run per file, concurrent
with the mutation battery on 10 cores, so upper bounds; "cold" is cold
for the first file only, because later files of the cold pass meet entries
earlier files wrote):

.. list-table::
   :header-rows: 1
   :widths: 52 16 16 16

   * - File (rows)
     - No memo
     - Memo, cold
     - Memo, warm
   * - ``test_trajectory_resolvent_reference.py`` (46)
     - 414.7 s
     - 424.3 s
     - 55.4 s
   * - ``test_l1_standoff_slab_cylinder.py`` (14)
     - 1949.5 s
     - 1920.0 s
     - 1063.0 s
   * - ``test_phase_c_crosscheck.py`` (9)
     - 1149.1 s
     - 511.5 s
     - 41.8 s
   * - ``test_unified_matvec_cylinder.py`` (32)
     - 870.8 s
     - 55.9 s
     - 53.6 s
   * - ``test_peierls_greens_function_mr.py`` (5)
     - 169.7 s
     - 161.1 s
     - 1.8 s
   * - ``test_peierls_greens_function_cylinder_mr.py`` (10)
     - 559.2 s
     - 566.2 s
     - 2.1 s
   * - ``test_reference_body.py`` (59)
     - 20.5 s
     - 14.5 s
     - 2.3 s

These are the prototype's numbers, before the review round and before the
exact medium left the clients; the specification's section 3 states the
measurement protocol the shipped memo owes (cold, warm and bypassed, three
runs each, serially). The cache after both passes held 8 724 KB.

**The gates** are 128 rows `[M]` (``--collect-only`` at ``bacf0787``): the manifest (22,
``test_traced_memo_manifest.py``), validation and the payload (18,
``test_traced_memo_validation.py``), the process and the store (27,
``test_traced_memo_process.py``), data files (1,
``test_traced_memo_data.py``), the clients (13,
``test_traced_memo_clients.py``), the three reviews' findings (33,
``test_traced_memo_findings.py``) and the ambient state (14,
``tests/gates/test_ambient_state.py``), all under ``tests/gates/``. The
real-tree witnesses (M4.4) copy every ``*.py`` under ``orpheus/`` into a
temporary directory and make a behaviour-neutral edit there: an edit to a
construction helper the reading ran makes the reading stale and leaves
the solve valid; an edit to the solver makes both stale, the reading
through its child pin; an edit to the cylinder's rays, which a sphere
reading never ran, leaves both valid; a constant in the reading module's
skeleton makes the reading stale.


Declared limits
===============

Each limit is a phase of #405 in the plan of record
(``.claude/plans/reference_cache.md``):

* **Three clients only.** No reference family besides the multi-region
  trajectory resolvent is memoised; the families migrate in phase P4, and
  a certification run becomes a client there. No gate yet refuses a
  traced memo on a withdrawn generator (the generating process refuses it,
  since its declared environment holds no ``ORPHEUS_*`` switch).
* **No cache on continuous integration.** The cache lives under
  ``.cache/references/`` on one machine; the continuous-integration job,
  the read-only test shards and the scheduled cold rebuild that would
  control the recorder itself are phase P5. A coarse cache key is safe
  there because a restored stale entry is a miss, never a wrong answer.
* **The shipped memo's suite timings are not measured.** The figures above
  are the prototype's; the protocol of the specification's section 3
  (cold, warm and bypassed, serially, three runs each) is owed.


What was tried and failed: the review round
===========================================

The memo landed in two commits (steps 1 to 3, then the clients) and was
reviewed by qa and the elegance-enforcer in parallel. qa reproduced 13
findings on the first commit, 7 of which served a stale or a wrong value;
the elegance-enforcer found 2 violations and 3 concerns, overlapping qa's
on the def walk, the type tree and the data pins. All were fixed in P3, at
their roots, before anything reached ``main``. Each finding is now a gate
in ``tests/gates/numerics/test_traced_memo_findings.py`` (18 rows), the
assertion form of qa's reproducer. The rows lettered A, P, E and N are
earlier: the test-architect's gates found them in the prototype the
specification was written against, before the memo was built. The three
unnumbered rows are not review findings: they arose in the fix round.

.. list-table:: The review's findings, their mechanisms and their fixes
   :header-rows: 1
   :widths: 6 40 34 20

   * - #
     - What failed, and how it hid
     - The fix, at the root
     - Gate
   * - qa 1
     - Same-name defs: the line found the def that ran, the pin stored the
       last def of that name. Stale hit, ``11.0`` against ``12.0``.
     - The pin is the qualified name and the def's rank.
     - ``test_q1_…``
   * - qa 2
     - Context variables are empty in a new thread: a memo called from a
       worker thread was not pinned and wrote to the default root. Stale
       hit, ``11.0`` against ``15.0``.
     - One recorder, process-wide; one process-wide target stack.
     - ``test_q2_…``
   * - qa 3
     - The payload codec admitted by ``isinstance``: ``NamedTuple``,
       ``IntEnum`` and masked arrays came back as their bases.
     - Exact types only; a subclass is refused.
     - ``test_q3_…``
   * - qa 4
     - The type tree stopped at the top level: ``sign(-0.0)`` served
       ``1.0``; ``Spec(8.0)`` served ``Spec(8)``'s entry.
     - The key is ``encode_exact``; the type tree retired.
     - ``test_q4_…``
   * - qa 5
     - Unrecorded dependencies: a missing file, a directory listing, a
       bytes path, a read relative to the working directory.
     - Data pins with ``"absent"``; listing pins; bytes decoded; the
       working-directory pin.
     - ``test_q5_…`` (2)
   * - qa 6
     - Defs under a module-level ``for``, ``while``, ``match`` or
       ``try``/``except*`` were never pinned. Stale hit, ``101.0`` against
       ``105.0``.
     - One walk descends every statement except a def's body.
     - ``test_q6_…``
   * - qa 7
     - A hit re-read the payload after verifying it; an interleaved write
       served a value no generation produced.
     - ``Hit`` carries the arrays it verified.
     - ``test_q7_…``
   * - qa 8
     - The memo writes after the recording stops, so entries written
       before a memo fix stayed hits.
     - ``MemoPin``.
     - ``test_q8_…``
   * - qa 9
     - Standard output was the result channel; a ``print`` in a generator
       made the caller fail.
     - A reserved descriptor; standard output redirected to standard
       error.
     - ``test_q9_…``
   * - qa 10
     - A directory entry replaced by ``rmtree`` and ``os.replace``: 458 of
       800 threaded writes failed and leaked staging directories.
     - One ``.npz`` file, one ``os.replace``.
     - ``test_q10_…``
   * - qa 11
     - (a) An exception that cannot be rebuilt arrived anonymous; (b) an
       undeclared dataclass return was written, then refused on every
       load.
     - (a) A described ``RuntimeError``; (b) refused at the write.
     - ``test_q11_…`` (2)
   * - qa 12
     - A process importing ``orpheus`` through the editable finder raised
       ``Unpinnable`` for every source file.
     - A file on no ``sys.path`` entry is pinned by its absolute path.
     - ``test_q12_…``
   * - qa 13
     - A self-call ``f(x)`` inside ``f(x)`` started interpreters without
       bound.
     - ``RecursionError`` on the ancestry.
     - ``test_q13_…``
   * - S1
     - The manifest coerced rows by position: a child pin placed among
       functions became a def pin (pins of equal arity cannot be told
       apart).
     - One set of tagged pins, parsed once in ``from_json``.
     - the manifest gates
   * - S2
     - Two tracers, already differing.
     - One ``Recording``.
     - the manifest gates
   * - S4, S5, C1 to C3
     - A second, shallower walk for the key; an untyped decorator and
       unchecked canonical keys; the constructor form spelled three times;
       two context variables.
     - ``encode_exact``; ``TracedMemo[P, R]`` with a ``ParamSpec`` and a
       checked ``canonical``; ``content.constructor_arguments``; one
       target stack.
     - —
   * - A
     - (prototype) ``typing.get_type_hints`` in the decoder ran a child's
       ``__annotate__``, reported at the child's def line, so the parent
       listed the child's def: M3.12b passed through the parent's own
       manifest, green for the wrong reason, caught by M3.12a's negative
       leg.
     - Signature scopes are pinned by the skeleton, never as a body.
     - M1.3c
   * - P
     - (prototype) ``<generic parameters of find_factor>`` has no def of
       its name at its line; refused as unpinnable, every generation
       failed.
     - The same rule.
     - M1.3c
   * - E
     - (prototype) The editable ``orpheus`` finder belongs to no
       top-level name; refused, every generation failed.
     - The ``RECORD`` search; an editable file pinned as source.
     - M1.5
   * - N
     - (prototype) ``isinstance(v, float)`` tested first wrote
       ``np.float64`` as ``float``.
     - Numpy scalars tested first, by exact type.
     - M2.7
   * - —
     - Directory listings made by the import system were pinned, so the
       M4.4 witnesses missed on every new file in a ``sys.path``
       directory.
     - Listings made by the import system are excluded.
     - M4.4
   * - —
     - Programs: ``platform.processor()`` runs ``uname`` during the clients'
       imports, an unrecorded dependency of every entry.
     - Started programs pinned by their bytes; a started Python refused.
     - ``test_a_started_program_…``
   * - —
     - The boot script ran as ``__main__`` and the memo imported a second
       copy with no active recordings (found while fixing).
     - The boot aliases itself in ``sys.modules``.
     - the process gates

**The batteries.** Before the fixes, qa re-targeted its mutation battery
on the real module: 37 of 38 arms reddened their target row (arm B14,
the non-atomic write, blind as the specification declared; arm A5
crash-dominated). After the fixes, a fix battery of 8 arms, one per
fixed mechanism (the def walk, the last-def pin, the exact key, standard
output, the thread children, the memo pin, the listings, the cycle guard),
reddened its gate in 8 of 8, each arm restoring the pristine files and
byte-comparing them after the run. Before any of this, the
test-architect's specification battery ran 42 arms against the
prototype, and 41 reddened their row; the blind one is the atomic write.

**Refuted by qa, kept as facts:** the ``.npz`` bytes are deterministic, so
regenerating a child with one answer does not make its parents stale;
the parent and the child agree on the ``-O`` level; a pin relative to a
``sys.path`` entry validates the file Python would import.


The second review
-----------------

qa reviewed the fixed memo (``39ee20f2``) again. 11 of the 13 first-round
findings were closed; two were half open (finding 5: an absence probed by
``exists()``; finding 8: ``content.py`` runs on the write path after the
recording stops); 8 findings were new, 6 of them serving a wrong value. The
spy rows of M4.6 had their teeth back (qa: "14 of 14 red under the
fixture; 15 of 15 green without it"), and qa's battery
reddened its target in 43 of 46 arms. Each finding was fixed at its root
in ``994ba740``, read against two standing rulings: the key separates what
the function can tell apart (Q2), and a stale value is never served. Ten
rows joined ``test_traced_memo_findings.py``.

.. list-table:: The second review's findings, their mechanisms and their fixes
   :header-rows: 1
   :widths: 8 40 34 18

   * - Finding
     - What failed, and how it hid
     - The fix, at the root
     - Gate
   * - r1
     - The key was exact at its leaves only: a masked array and its plain
       twin, two dictionary orders, a sparse matrix's explicit zero and
       its format, a ``NamedTuple`` and its tuple, and a ``compare=False``
       constructor field each shared one key and served the other call's
       answer.
     - The key encodes what the generating process receives: each object
       by its constructor form, each container by its concrete type, each
       mapping in order, each sparse matrix by its storage, each subclass
       of ``str``, ``bytes`` or ``ndarray`` by its type and state.
     - ``test_r1_…``
   * - r2, r3
     - A spawn-context process pool ran Python no recording saw (its start
       raises ``_posixsubprocess.fork_exec``, not ``subprocess.Popen``);
       a started program was pinned by its launcher's bytes while the
       files it read went unpinned: ``cat`` served ``1.0`` against
       ``50.0``, and ``shell=True`` and ``env python3`` slipped past the
       refusal of a started Python.
     - Only a declared machine query (``uname``) may be started; every
       other start is refused, ``fork_exec`` included.
     - ``test_r2_r3_…``
   * - r4
     - A directory listed by a module body while it was being imported was
       exempt as the import system's (an import was on the stack): ``11.0``
       served against ``22.0``.
     - The exemption is decided by the immediate caller.
     - ``test_r4_…``
   * - r5
     - Every ``.py`` read was skipped as code, so a generator reading a
       ``.py`` file as data served ``3.0`` against ``30.0``.
     - The blanket skip retired; only the import system's own reads are
       exempt.
     - ``test_r5_…``
   * - (environment)
     - A generator reading a variable outside every key served ``2.0``
       against ``3.0``.
     - A declared environment, pinned; ``PYTHONHASHSEED`` fixed. The
       ``ORPHEUS_*`` scrub retired, subsumed.
     - ``test_m3_16_…``
   * - r7
     - The working-directory pin held the run's FINAL directory, so a
       generator that changed directory and read a relative path was
       stale on every call (its own fresh entry did not validate).
     - The pin holds the directory the run started in.
     - ``test_r7_…``
   * - (after the round)
     - A parent pinned its child's source files as DATA: the lookup of the
       child inside the parent's generation validated the child's manifest
       by reading its sources, and those reads were recorded. `[M]`
       2026-10-04 at ``994ba740``: the reading entry's 105 ``.py`` data
       pins were exactly the solve entry's 105 module pins (105 of 105 in
       both directions), and 368 presence pins against 55 after the fix; a
       comment edit to any of the solve's modules made every reading stale
       (a spurious miss, never a stale hit). The M4.4 witnesses missed it:
       the file their "other geometry" row edits is not among the solve's
       modules.
     - The memo's lookups run quiet (``unrecorded``); the parent holds its
       child by its child pin only (``6bfde984``).
     - ``test_a_parent_does_not_pin_its_childs_sources_as_data``
   * - (threads)
     - A ``bypass()`` in one thread redirects a concurrent memo call in
       another.
     - Declared, not fixed.
     - —
   * - q5 (half)
     - An answer that depended on ``exists()`` finding no file was served
       stale: ``1.0`` against ``99.0``.
     - ``os.stat`` and ``os.lstat`` probed while recording;
       ``PresencePin``.
     - ``test_q5_absence_probed_by_exists_…``
   * - q8 (half)
     - ``content.constructor_arguments`` runs on the write path, after the
       recording stops, and was not pinned.
     - ``MemoPin`` covers ``content.py``.
     - ``test_q8_the_memo_pin_covers_the_write_path``
   * - q10 (reader)
     - The writer-only gate was blind to an in-place write.
     - A reader leg: while four threads rewrite one key, every lookup is
       ``Hit`` or ``Absent``, never ``Corrupt``.
     - ``test_q10_a_reader_…``
   * - q13 (bound)
     - With the cycle guard removed, a chain through distinct keys hung
       instead of failing.
     - Nested generations bounded at 16.
     - ``test_q13_a_chain_…``

**Found while fixing, not by the review:**

* **The audit event** ``exec`` **is the builtin** ``exec()``, which builds a
  dataclass's generated methods; it is not a process start, so it is not
  among the spawn events.
* **The exact key refused** ``EmissionSpectrum``, the stateless
  ``ndarray`` subclass held by ``Mixture.chi``, once the key refused every
  ``ndarray`` subclass: `[M]` 23 rows of the affected suites were red. A
  subclass is now its type and its own instance state beside its bytes, so
  ``EmissionSpectrum`` is its type and its bytes, and only a subclass whose
  state has no content identity (a masked array) is refused.
* ``importlib.metadata`` **lists every** ``sys.path`` **directory** to find
  the installed distributions, and with the exemption limited to the frozen
  bootstrap those listings pinned the working directory. The import system
  is now the frozen bootstrap AND the ``importlib`` package.
* **The import-system check recursed.** It ran ``realpath`` inside the
  ``os.stat`` probe, and ``realpath`` calls ``os.stat``. The check now
  compares raw file names against prefixes computed once, at import.

`[M]` 2026-10-04 after the fixes: the affected suites (numerics, reference,
data, mesh, geometry, the trajectory and multi-region solver tests), 7 487
passed, 0 failed.

**An instrument lesson.** qa's own reaper for orphaned generating
processes, ``pkill -f _traced_memo_boot``, matched by process NAME, so it
also killed the generating processes of a full-suite run going on at the
same time: all 5 failures of that run at ``39ee20f2`` were the reaper's,
and none was the memo's. A battery's orphans are killed by process tree,
never by name. `[M]` 2026-10-04, the full non-slow suite at ``994ba740`` in
a quiet worktree: 15 223 passed, 0 failed.


The third review, and the contract
----------------------------------

qa reviewed ``994ba740`` and ``6bfde984`` a third time and found further
ways an ADVERSARIAL generator can hide a dependency: a file's size or
link-ness read through ``os.stat``, where the presence pin keeps only the
file's kind; a data read through ``pkgutil.get_data``, behind
import-system frames; an optional import's absence; nine exotic argument
types that shared a key (K1 to K9); ``uname`` admitted by its basename;
and ``PATH`` in the environment pin. None of them touches the three real
clients. The user ruled the close (the section on the generator contract
above): the cheap, real fixes land in ``bacf0787``, and the rest is
declared.

.. list-table:: The third review's findings and what was done
   :header-rows: 1
   :widths: 30 50 20

   * - Finding
     - Fixed, or declared
     - Gate
   * - The exact key (K1 to K9): set iteration order; a sparse matrix's
       class beside its format; a sparse matrix with no stored arrays; a
       ``defaultdict``'s factory; an ``int`` subclass whose ``__str__``
       lies; ``__slots__`` state; an ``InitVar``
     - Fixed: set order kept; sparse class keyed and a ``dok_matrix``
       refused; mappings other than ``dict`` and ``FrozenMapping``
       refused; ``int`` and ``float`` subclasses keyed by type and state,
       the value never through ``__str__``; ``__slots__`` in subclass
       state; ``constructor_arguments`` refuses an ``InitVar`` class,
       which also refuses it at the process boundary and in the payload.
     - ``test_k_…``
   * - ``uname`` admitted by its basename: a script named ``uname``
       earlier on ``PATH`` was admitted
     - Fixed: admitted by its system path; ``fork_exec`` names the program
       by its executable candidates.
     - ``test_u_…``
   * - ``PATH`` pinned: activating a virtual environment made every entry
       stale; and a generator that set a declared variable never validated
     - Fixed: the pin holds the declared variables, ``PATH`` excluded, as
       the run received them.
     - ``test_e_…``
   * - File size or link-ness through ``os.stat``; a read through
       ``pkgutil.get_data``; an optional import's absence
     - Declared: outside the generator contract; P5's cold rebuild is the
       witness.
     - —

Development history
===================

Reverse-chronological changelog of phase P3 of #405, the plan
``.claude/plans/reference_cache.md`` and its gate specification
``.claude/plans/reference_p3_spec.md``. Every row below sits on the branch
``refactor/p3-ambient-state`` at the time of writing;
``git merge-base --is-ancestor <hash> main`` outranks this column.

.. list-table::
   :header-rows: 1
   :widths: 11 55 10 24

   * - When
     - Milestone
     - Issue
     - Where
   * - 2026-10-04
     - **P3 closed by the generator contract** (the user: "Contract +
       cold-rebuild control"), a scope boundary whose witness is P5's cold
       rebuild; the third review's cheap fixes (the exact key's set order,
       sparse class, mapping, ``int``/``float`` subclass, ``__slots__``
       and ``InitVar`` cases; ``uname`` by its system path; the
       environment pin without ``PATH``, as received). Three gates (K, E,
       U).
     - #405
     - ``bacf0787``
   * - 2026-10-04
     - **This chapter after the second review**, and the regenerated
       matrix.
     - #405
     - ``ab2dcd66``, ``889cae2f``
   * - 2026-10-04
     - **The memo's own store reads are not a generation's dependencies**:
       a parent no longer pins its child's sources as data (found by this
       chapter's manifest measurement).
     - #405
     - ``6bfde984``
   * - 2026-10-04
     - **The second review round.** The exact key exact all the way up;
       only ``uname`` may be started; the import-system exemption by the
       immediate caller, over the frozen bootstrap and ``importlib``;
       presence probes (``PresencePin``); a declared, pinned environment
       (``EnvironmentPin``); the starting working directory; ``MemoPin``
       over ``content.py``; nested generations bounded at 16. 10 finding
       gates; the affected suites 7 487 passed.
     - #405
     - ``994ba740``
   * - 2026-10-04
     - **This chapter**, and the regenerated verification matrix.
     - #405
     - ``51c8485b``, ``bfc926eb``
   * - 2026-10-04
     - **The review round: no stale or wrong value served.** qa's 13
       findings (7 serving a stale or wrong value) and the elegance
       review's, each fixed at its root: one def walk with ranks, nine
       tagged pin kinds, one process-wide recorder, ``encode_exact`` as
       the key, exact payload types, one atomic file, a reserved result
       descriptor, the cycle guard, the boot module aliased. 18 finding
       gates; a fix battery reddened 8 of 8.
     - #405
     - ``39ee20f2``
   * - 2026-10-04
     - **Step 4: the clients.** The two multi-region solvers memoised with
       their six array parameters canonical;
       ``TrajectoryResolventDerivation`` gains content identity and its
       ``evaluate`` is a memoised method; the exact infinite medium is not
       a client (Q5); the solve-spy rows read under the bypass.
     - #405
     - ``330ba98a``
   * - 2026-10-04
     - **Steps 1 to 3: the traced memo**, written from the
       test-architect's prototype: the manifest, validation, the payload,
       the key, the process boundary, the store and the bypass.
     - #405
     - ``99897d10``
   * - 2026-10-04
     - **Step 5: ambient state.** The slab route is an argument
       (``SlabRoute``), and the two global ``mp.dps`` writes are local
       ``workdps`` contexts.
     - #405
     - ``a0f1b6ef``
   * - 2026-10-04
     - **The specification ruled**: 5 steps, 96 gate rows, a 42-arm
       battery on a prototype; the user's rulings Q1 to Q9, Q5 excluding
       the exact medium.
     - #405
     - plan ``c962dada``
   * - 2026-10-04
     - **P3 opened**: the measurements, the two architectures, and the
       user's four rulings (architecture B, the reading and its solve
       child, the three clients, docstrings stripped); the census.
     - #405
     - plan
