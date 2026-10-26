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
   nothing; what a validation keeps for the process; the source index and
   its write rule; what a miss costs; the
   exact payload; the process boundary; the bypass; how the
   memo relates to the reference certificate; the exact medium's
   exclusion; the clients; the limits no manifest can see; the review
   round's defects and their fixes; the ambient-state step (the slab
   route as an argument, the mpmath precision as a local context).
   Not owned here: the ontology of a reading, a claim and a certificate
   (readings_and_certificates.rst, verification-reference-architecture);
   the content encoder's canonical forms
   (structured-geometry-content-identity); the trajectory resolvent's
   reading (references/characteristic_origins.rst,
   characteristic-origins-history-reading); the slab routing's physics
   (references/peierls_nystrom.rst, theory-peierls-slab-polar-g5-routing).
   Plan of record: .claude/plans/reference_cache.md, from "P3 opened
   (2026-10-04)" to its end, and discussion 2; specification of the
   gates: .claude/plans/reference_p3_spec.md. Written at P3's close,
   2026-10-04, against branch refactor/p3-ambient-state at 39ee20f2;
   updated for the second review round at 994ba740, for the generator
   contract and the third review at bacf0787, and for #602 (the miss
   cost, the pinned-source loader, per-file ownership, the source index)
   at 8e02236f, e666e58d, 0c569df8 and 155fbdbd, and for cut E (the fork
   template, every input observed at its first touch, the RECORD pin) on branch
   perf/traced-memo-fork-template, 2026-10-16, and for the hit cost (a
   validation's reads kept for the process: FileIdentity, Witness, Kept)
   at 8fbb1d61 on branch fix/memo-hit-cost, 2026-10-17, and for the
   first test-oracle clients and the memo's child resolution
   (_resolve_memo) at 7ca91c43 and 333268fb on branch
   test/characteristic-oracle-memo, 2026-10-17, and for the two oracle
   routes stored after them at 76d0409d on branch test/oracle-memo-columns,
   2026-10-17, and for the slow cylinder point readings read through the
   door (concept 2) at 94143ecd on branch test/door-cylinder-readings,
   2026-10-17, and for the PS-1982 vacuum-sphere solver as the eighth
   client (concept 3) on branch test/ps1982-memo over d5c1ea54,
   2026-10-17.

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
interpreter (or, where forking is measured safe, in a child forked from a
template interpreter that holds nothing a fresh one would not) that records every function that starts, every file it opens
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
  generation per miss: the interpreter start and the imports, about
  0.5 s, plus the computation and the manifest, which bring a cheap miss
  to 0.66 to 0.85 s against a 0.010 s hit (`[M]` 2026-10-15 at
  ``155fbdbd``, five characteristic transport blocks that take 0.02 to
  0.17 s in process, three repeats, a warm bytecode cache and warm source
  indexes). Before #602 the same miss cost 2.65 to 2.88 s; the section
  "What a miss costs" below attributes it.
* **A fresh process compiles from the bytes it pins.** A loader the boot
  script installs before any ``orpheus`` import compiles every
  first-party module from the bytes it reads, and the manifest pins
  exactly those bytes, never the file read again after the run. A
  timestamp ``.pyc`` validated by a whole-second modification time and a
  size had let a regenerated entry hold the old answer and validate
  (ERR-115; the section "A fresh process must compile from the bytes it
  pins").
* **Each source text is parsed once per content.** What a manifest and
  a validation read from a first-party source file (its skeleton digest,
  every outermost def's span and digest, the span of every node a code
  object is compiled from) is kept under ``.cache/sources/``, named by the
  text's bytes, the interpreter and the memo's own source, and only a
  generating process, which knows the bytes of the memo it ran, writes
  one (the section "The source index" below).
* **The key is exact, finer than content equality.** It is
  :func:`~orpheus.values.content.encode_exact` of the function's
  identity, its signature-bound arguments (defaults applied, each declared
  canonical form put in) and the platform tag. ``8`` and ``8.0``, ``-0.0``
  and ``0.0``, an integer array and its float twin are different keys,
  because a function can tell them apart even though ``==`` cannot; and
  the key is exact all the way up: an object is its constructor form, a
  container carries its concrete type, a mapping its order.
* **Validation runs nothing.** A lookup re-hashes the manifest's pins
  against the current files, reading again only a file or directory whose
  identity on disk (device, inode, size, modification and change times)
  has moved since the process last read it, and returns one of four verdicts,
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
* **Eight clients today.** Three are in the characteristic reference: its
  reading, ``CharacteristicDerivation.evaluate``, keyed on the derivation
  and the observable; its solve, ``CharacteristicDerivation.solve``, keyed
  on the derivation alone and a child entry of every reading (#592); and
  one group's transport block, ``group_transport``, keyed on the basis,
  the walls, the group's total cross sections, the resolution and the
  emitting regions, and a child entry of the solve (#405 P3 step 2c). Four
  are a test's oracle, the mpmath routes ``phi_cylinder_vacuum``,
  ``phi_cylinder``, ``specular_cylinder_total`` and
  ``cylinder_polar_density`` of
  ``tests/gates/derivations/_characteristic_point_mp.py``, which take
  seconds to minutes a call (#405 P3 concept 1 and ``76d0409d``,
  2026-10-17). The slow cylinder point-reading rows of the characteristic
  gates read through the first of these, the reading, since concept 2
  (``94143ecd``), which added no client. One is the independent reference
  of the characteristic reference's vacuum sphere, the PS-1982 solver
  ``solve_ps1982_vacuum_sphere``, 11.6 to 20.7 s a solve (#405 P3 concept
  3, 2026-10-17). P3 shipped three others, the two multi-region trajectory-resolvent solvers
  (``solve_greens_function_sphere_mr`` and
  ``solve_greens_function_cylinder_mr``) and the trajectory-resolvent
  derivation's ``evaluate``; they were deleted with their family at P1
  step (e2) of the characteristic reference campaign (2026-10-10). The
  exact infinite medium is **not** a client: its reading costs 0.01 s,
  less than one interpreter start (the user's ruling, Q5).
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

A fresh process removes what memory could hide, and the section
:ref:`verification-reference-cache-pinned-source-loader` below states
what it does not remove on its own: the code it runs must be compiled
from the bytes its manifest pins.

Since #602's cut E a miss may be generated by a child forked from a
**template** interpreter instead of a new one, where forking has been
measured safe. The ruling still holds, because the template is built to
hold nothing a fresh process would not: it is itself a fresh interpreter,
never a fork of the asking process (hazard 3), and it runs no generation,
only the imports, so no generation inherits another's memory (hazards 1
and 2). How, and the gates that pin it:
:ref:`verification-reference-cache-fork-template`.

**What crosses into the process is content, never a built object.** The
census of 2026-10-04 (finding F1) measured the trap: unpickling a built
``TrajectoryResolventDerivation`` or ``Billiard`` (both deleted since)
restored the object
without running its construction, so in 2 of 2 traced trajectory runs the
construction functions (``Billiard.__post_init__``, ``_route``,
``_layered_xs_payload``, ``_read_isotropically``, ``reference_body`` and
the derivation's ``__post_init__``) were absent from the record; in the
positive control, which constructs the derivation inside the process, 7 of
7 were present. So every dataclass crosses through its constructor
(``_ConstructorPickler`` rebuilds it from
:func:`~orpheus.values.content.constructor_arguments`, and a
:class:`~orpheus.values.content.ContentIdentity` value already pickles
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


.. _verification-reference-cache-pinned-source-loader:

A fresh process must compile from the bytes it pins
---------------------------------------------------

A fresh process has nothing in memory, but it does not follow that it
runs the source files as they stand on disk. CPython does not compile an
imported module from its source when it finds a valid ``.pyc`` in the
module's ``__pycache__``, and the default ``.pyc`` (the timestamp form of
PEP 552) is valid when two fields of its header match the source file:
the source's modification time truncated to whole seconds, and its size.
A same-size edit whose modification time falls in the same whole second
as the write the ``.pyc`` recorded passes that check, and the interpreter
runs the bytecode of the text before the edit.

Until ``e666e58d`` the manifest pinned a first-party source file by
reading it after the run. The two definitions of "the source of the code
that ran" then disagreed, and the memo served a wrong answer as a hit:
the parent's import wrote the ``.pyc``, a same-size edit with its
modification time in the same second made the next lookup ``Stale``, the
generating process ran the old bytecode, and the manifest pinned the new
text. `[M]` 2026-10-15, the test-architect's probe
``scratch/characteristic_architecture/p3/m602/ta/p8_pyc_through_memo.py``:
a generator returning ``x * SCALE``, edited from ``SCALE = 3.0`` to
``SCALE = 4.0``, regenerated **4.5** for ``generate(1.5)`` where **6.0**
is right, and the next lookup was a ``Hit``. The same re-read was wrong
without any ``.pyc``: an edit landing during a generation, after the
module was compiled and before the manifest was written, was pinned
beside an answer computed from the text before it. The defect is
ERR-115 in :ref:`theory-verification-error-catalog`, with how it hid.

**The repair is a loader, not a check.** The generating process compiles
every first-party module from bytes it reads itself, and the manifest
pins exactly those bytes (the orchestrator's ruling, 2026-10-15). The
manifest then cannot pair an answer with a text other than the one that
produced it, whatever the file's times are and whenever an edit lands.
Comparing the ``.pyc`` against the source a second time would have
narrowed the window without closing it, since the edit during a
generation needs no ``.pyc``.

The parts, all in :mod:`orpheus.numerics._traced_memo_boot`:

* **What it covers.** The boot script starts the recording on its first
  line and, before any ``orpheus`` import, calls
  ``_install_pinned_sources``. From then on every module whose spec
  carries a plain ``SourceFileLoader`` and whose origin lies outside the
  standard library and the site-packages directories is compiled by
  ``_PinnedSourceLoader``. Code in the standard library is pinned by the
  interpreter and code in site-packages by its distribution's version,
  never by its bytes, so it keeps the ordinary loader. An extension
  module, a namespace package and a sourceless module are left as they
  are (only a spec whose loader is exactly ``SourceFileLoader`` is
  swapped).
* **The finder is first and delegates.** ``_PinnedSourceFinder`` is
  inserted at ``sys.meta_path[0]``. It asks every other finder on
  ``sys.meta_path``, in order, for the module's spec, and swaps the
  loader of the spec the first one returns. A path hook alone would see
  only modules found on ``sys.path``, and an editable install resolves
  ``orpheus`` through its own meta-path finder when ``sys.path`` does not
  hold the package; a finder appended after that one would never be
  asked. (Gates P4, a module only a meta-path finder knows, and P4b, a
  process importing ``orpheus`` through the editable finder. `[M]` the
  battery arm that appends the finder instead of inserting it reddens P4
  among 58 rows.)
* **The bytes it records.** ``get_code`` reads the source once, offers
  the bytes to every active recording whose ``sources`` field is a
  dictionary (keyed by the file's real path), and compiles from them.
  ``_compiled_bytes`` in :mod:`orpheus.numerics.traced_memo` gives the
  manifest those bytes, so the skeleton and every def digest are of the
  text that ran.
* **Its reads are the memo's, not the run's.** The whole body runs
  under ``unrecorded()``. Recorded, its reads of the source would pin
  every imported module as data, so a comment edit would make the entry
  stale, and its read of the bytecode cache would pin a cache file
  (gate P3; `[M]` the arm that drops ``unrecorded()`` reddens P3).
* **The bytecode cache** lives at
  ``.cache/bytecode/<cache_tag>-O<level>/`` under the repository root,
  beside the reference store and the source indexes (one root, ``CACHE``
  in the boot module), one directory per interpreter and ``-O`` level. An entry is named by the hexadecimal BLAKE2b-256 digest of the
  module's path, the interpreter's cache tag, the ``-O`` level and the
  source bytes. Its content is that key (32 bytes), the BLAKE2b-256
  digest of the payload (32 bytes) and the payload, the code object
  serialised by ``marshal``. The ``-O`` level is part of both the
  directory and the key because ``-O`` acts at compile time: it makes
  ``__debug__`` a false constant the compiler folds and removes every bare
  ``assert``, so one source compiles to two different code objects, and
  an entry shared between the levels would serve one level's code to the
  other (gate P6, a generator returning
  ``x + (1.0 if __debug__ else 2.0)`` generated at both levels; `[M]` the
  arm that fixes the level at 1 in the directory and the key reddens P6).
* **A corrupt or foreign entry is recompiled, never run.** An entry is
  run only if its first 32 bytes equal the key computed from the bytes
  just read and its second 32 bytes equal the digest of the rest. Any
  read error, a short file, a flipped bit in the payload or in the stored
  digest, and another file's entry copied under this name all fall
  through to a compile (gate P9, five rows; `[M]` the arm that replaces
  both comparisons with a length check reddens the stored-digest row and
  the copied-entry row). The key carries the source bytes, so an entry of
  an old text is never found for a new one: the cache cannot be stale,
  only corrupt, and corruption is detected.
* **The write is atomic.** A compiled entry is written by
  ``write_atomically``, the one atomic write of the memo's three stores,
  to ``<entry>.<pid>.<thread id>.part`` and moved into place with
  ``os.replace``, so a concurrent reader finds a whole entry or none,
  never a part (gate P8: two processes each delete the entry and compile
  one module 200 times while the test process reads the entry as fast as
  it can, and every read that finds the file finds a whole entry; the
  file is the instrument, because the loader recompiles a partial entry
  silently and no answer can show it; `[M]` the arm that writes straight
  to the entry's name reddens P8). A write that fails (a read-only checkout) is ignored: the
  module was compiled, and the next process compiles it again.
* **Why apart from** ``__pycache__``. The ``__pycache__`` files are the
  ones validated by timestamp, which is the hazard; they are also what
  every ordinary import, the test process's included, reads and writes,
  so a generation leaves them alone (gate P5: a module only the
  generation imports leaves no ``.pyc`` in its ``__pycache__``, and its
  code is in the bytecode cache under its key, in the format above). The
  entry is not a ``.pyc`` of any form: the round-two specification had
  named CPython's checked-hash ``.pyc``, and the shipped format carries
  its key and the digest of its payload instead.

**Code the loader did not compile is refused.** In a generating process a
code object whose file is a first-party source the loader never compiled
(an ``exec`` of a file read by hand, ``runpy``) has no recorded bytes, and
a re-read after the run is the very pin this section removes, so
``_compiled_bytes`` raises
:class:`~orpheus.numerics.traced_memo.Unpinnable`, naming exec and runpy
(gate P7). `[M]` the test-architect's AST census over the 390 files of
``orpheus/`` (``ta/p9_exec_census.py``, its positive controls found all 8
spellings): 4 calls of ``exec``, ``eval``, ``compile``, ``runpy.*``,
``spec_from_file_location``, ``exec_module`` or ``SourceFileLoader``, and
0 of the 4 compile a first-party file (two are ``re.compile``, and
``numerics/mesh_free_function.py`` evaluates a syntax tree under the
filename ``<srepr>``, which the manifest drops). No generator relies on
the refused case.

**What it costs.** `[M]` 2026-10-15
(``scratch/characteristic_architecture/p3/m602/run_after_loader.log``
against ``run_after_abcd1.log``, the fixture of
:ref:`verification-reference-cache-miss-cost`, three repeats, a warm
bytecode cache): the median time to spawn and finish a generation rose
by 3 to 35 ms per block, 21 ms on average over the five blocks. Reading
the source and computing one BLAKE2b digest per module is the added
work; the compile itself is skipped by the bytecode cache as the
``.pyc`` skipped it before. The cost of a cold bytecode cache, a
compile of every first-party module the generation imports, is not
measured.


.. _verification-reference-cache-key:

The key: exact, not content equality
====================================

The key is the 256-bit BLAKE2b digest of

.. code-block:: text

   encode_exact((function id, signature-bound arguments, platform tag))

where the function id is ``module:qualname``.

**The arguments are bound to the signature**, defaults applied
(``inspect.signature(f).bind(...).apply_defaults()``). Without the binding
one call has two keys: ``Billiard`` (deleted since) omitted the solver
settings it left at ``None``, while a direct caller may pass the same
default explicitly (census, 2026-10-04).

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
Q2). The content encoder :func:`~orpheus.values.content.encode` follows
``==`` by the user's ruling of 2026-10-02
(:ref:`structured-geometry-content-identity`): two values that are the
same physics are one value. A function can still tell such values apart,
so a key that followed ``==`` could serve one call another call's answer.
:func:`~orpheus.values.content.encode_exact` is the same walk of the
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
:mod:`orpheus.values.content`), rather than a second, shallower
description of the arguments (review finding qa 4). The exact policy
decides more than the leaves:

* **an object is its constructor form**,
  :func:`~orpheus.values.content.constructor_arguments` (or the
  arguments of its ``__reduce__``), which is exactly what the generating
  process receives, since every argument crosses through its
  constructor. A field declared ``compare=False`` that the constructor
  takes is therefore in the key, although content equality ignores it;
* **a container carries its concrete type**: a ``NamedTuple`` is not the
  plain tuple of its values;
* **a mapping keeps its order**:
  :class:`~orpheus.values.content.FrozenMapping` documents its iteration
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
* **a set is ordered by its elements' encodings**, as content identity
  orders it, and not by its iteration order (#592). A set's iteration
  order follows the process's hash seed: ``frozenset({1, 9})`` iterates
  ``[1, 9]`` and ``frozenset({9, 1})`` iterates ``[9, 1]`` in one process,
  and a frozenset of strings iterates differently under different
  ``PYTHONHASHSEED`` values. A key that kept iteration order therefore
  keyed one call differently in two processes. Every generation runs
  under ``PYTHONHASHSEED=0`` (:ref:`verification-reference-cache-boundary`),
  while the calling process (pytest, a user's session) runs under a
  random seed, so a parent and the generations it starts disagreed on the
  key of one call. `[M]` 2026-10-10, qa's #592 review: the solve key of
  one characteristic derivation, an eigen question on a body with an
  (n,2n) material, whose gauge ``CellCoefficient`` holds a fission cell
  and an (n,2n) cell in a frozenset, took two values over five unseeded
  processes (three and two); under a non-zero parent seed the parent's
  ``answer`` missed the solve entry its readings had written and wrote a
  second one. The earlier argument for keeping the order, that a function
  iterating a set can tell two equal sets apart, does not hold across the
  process boundary: the generation never receives the caller's order, it
  rebuilds its arguments under seed 0. Iteration order is therefore no
  input of a memoised function, and the generator contract says so
  (:ref:`verification-reference-cache-contract`). The gates are
  ``test_s1_one_call_is_one_key_under_every_hash_seed``, which keys a
  64-string frozenset and that derivation's solve in fresh interpreters
  under seven hash seeds and requires one key each (its positive control:
  the 64 strings iterate in at least two orders over the seeds), and its
  witness ``test_s1_the_iteration_order_arm_splits_the_keys``, which
  restores the iteration-order spelling in each child and requires the
  keys to split (``tests/gates/numerics/test_traced_memo_findings.py``);
* **a sparse matrix is keyed by its class as well as its format** (a
  ``csr_matrix`` is not a ``csr_array``), and one that stores no arrays (a
  ``dok_matrix``) is refused;
* **a mapping other than** ``dict`` **and**
  :class:`~orpheus.values.content.FrozenMapping` **is refused**: a
  ``defaultdict`` holds a factory beside its items;
* **a dataclass with an** ``InitVar`` **is refused**: its constructor
  takes a value it does not keep, so it has no constructor form.
  :func:`~orpheus.values.content.constructor_arguments` refuses it, which
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
       by its installed version and the digest of its ``RECORD`` text.
     - A dependency upgrade (numpy, scipy, mpmath, SymPy), and a rebuild
       reinstalled under the same version.
     - The design; the ``RECORD`` from qa's review of #602 cut E
       (finding 6).
   * - :class:`~orpheus.numerics.traced_memo.InterpreterPin`
     - The interpreter's version and bytecode tag.
     - The standard library and generated code changing under a new
       Python.
     - The design.
   * - :class:`~orpheus.numerics.traced_memo.DataPin`
     - A file the run opened for reading, by absolute path and the bytes
       it held when the run first opened it; ``"absent"`` when no regular
       file was there; a ``.py`` file read as data included. Also an admitted program the run started (``uname``), by
       its executable's bytes.
     - An edited data file; a file appearing where the run found none; a
       replaced program.
     - Q4 (M3.17); absence from qa finding 5; ``.py`` as data from the
       second review's finding r5; programs from the fix rounds; the
       first-touch observation from ERR-116.
   * - :class:`~orpheus.numerics.traced_memo.PresencePin`
     - A path the run probed without opening it (``os.stat``,
       ``os.lstat``, and through them ``os.path.exists``, ``isfile``,
       ``isdir`` and ``pathlib``'s ``exists``), by what was there at its
       first probe: ``"file"``, ``"directory"``, ``"other"`` or
       ``"absent"``.
     - An answer that branched on whether a file exists.
     - qa finding 5, the half the first fix left open.
   * - :class:`~orpheus.numerics.traced_memo.ListingPin`
     - A directory the run listed, outside the import system, by the
       sorted names in it at its first listing (or ``"absent"``).
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
       constructor form the write path runs; in a generation, by the bytes
       the pinned-source loader compiled them from.
     - An entry written by a memo whose writing code has since been
       fixed.
     - qa finding 8, both halves; the compiled bytes from the elegance
       review of #602 round two (F4).

**The normalised source of a def.** A def is hashed as the dump of its
syntax tree with line and column attributes excluded and every docstring
removed, so a comment, a layout change or a docstring edit never moves the
digest (the user's ruling of 2026-10-04, question 4: a docstring edit
cannot change an answer). Hashing source rather than bytecode lets
validation parse the current files without importing them. The source a
generating process hashes is the bytes it compiled, recorded by its
loader, never the file read again after the run
(:ref:`verification-reference-cache-pinned-source-loader`).

**A body that was only a docstring reads as** ``pass``, placed at the
docstring's lines, so the first line of such a def's body is still a line
of the file. ``_def_pin`` reads that line to tell a code object generated
from the def's signature (its ``__annotate__``) from one of its body.
Until ``e666e58d`` the substituted ``pass`` carried no location, and
``_def_pin`` raised ``AttributeError`` on the ``__annotate__`` code object
of such a def, the case of a ``typing.Protocol`` method stub whose
annotations a generator reads (`[M]` 2026-10-15, the test-architect's
``scratch/characteristic_architecture/p3/m602/ta/p1_pyc_and_pass.py``;
present since ``99897d10``). The generation failed loudly; no entry was
written. Gate A1 is that generation, which must succeed with a
``ModulePin`` and no ``DefPin`` for the stub, whose body never ran. Gate
A2 shows the location moves no digest: over every tracked Python file
with a docstring-only body, the strip's digests equal those of the same
text with each such docstring replaced by ``pass``. `[M]` 2026-10-15,
``scratch/characteristic_architecture/p3/m602/a2_arm.py``, in process: A2
is green on the tree, and red when the strip substitutes
``Expr(Constant(None))`` in place of ``pass`` (38 of 38 such files digest
unlike ``pass``).

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
the version of the distribution that OWNS it, the one whose ``RECORD``
lists the file; everything else is source, pinned by content. A file of an
**editable** distribution is pinned as source, because an editable
install's version never moves: every generating process runs the editable
``orpheus`` finder, ``__editable___orpheus_0_1_0_finder.py``, which the
``RECORD`` of ``orpheus`` lists. A ``site-packages`` file that no
``RECORD`` lists is refused as unpinnable; so is every file of an
``*.egg-info`` install, whose file list is not a ``RECORD`` (`[M]`
2026-10-15: no such install in the environment). :func:`~orpheus.numerics.traced_memo.origin`
takes the ownership table as a required argument, so no caller can build
it once per file.

Ownership is a property of the file, never of its top-level name. Until
2026-10-15 (#602) the owner was found through
``importlib.metadata.packages_distributions``, a map from each top-level
name to the distributions that claim it, and the first claimant was taken.
A namespace top is claimed by several distributions, so
``sphinxcontrib/tikz.py`` was pinned to ``sphinxcontrib-applehelp``, and an
upgrade of ``sphinxcontrib-tikz`` left valid every entry that ran its code
(`[M]` 2026-10-15, the test-architect: four namespace tops in the
environment, ``sphinxcontrib``, ``jaraco``, ``opentelemetry`` and
``_binding``; no client ran their code). The table, ``_Installed``, reads
every distribution's ``RECORD`` once per checkout (one manifest or one
validation, which holds it as ``_Checkout.installed``), as CSV because a
path may be quoted and hold a comma: 30 446 paths in 18 to 20 ms, where
the top-level map cost 0.22 s. Ownership stays per checkout. What is
installed as each distribution (``_Installed.installed``: its version and
the digest of its ``RECORD``) is answered from a read kept for the whole
process, ``_DISTRIBUTIONS``, made again when a ``sys.path`` directory or
the found distribution's ``METADATA``, ``PKG-INFO`` or ``RECORD`` changes
on disk (since 2026-10-17;
:ref:`verification-reference-cache-kept-reads`). So a checkout is two
parts: its first-party texts (``_Checkout``) and its installed
environment (``_Installed``: which distribution owns a file, and which
install of it is there).

**A distribution is pinned by its version and its** ``RECORD``. A
version alone does not name an install: a distribution rebuilt and
reinstalled under the same version (a local build, a patched wheel)
installs other files, and an entry that ran its code stayed valid. The
``RECORD`` lists the hash of every file the distribution installed, so its
digest moves when any of them does. qa's review of #602 cut E (finding 6)
showed it on a fake site directory; gate
``test_qa6_a_distribution_rebuilt_under_the_same_version_is_another_install``
rewrites the ``RECORD``'s hashes under the same version and requires the
manifest to read stale; battery F arm R-6, the comparison on the version
alone, reddens it alone (`[M]` 2026-10-16). The digest is of the ``RECORD``
TEXT as ``importlib.metadata`` reads it, in universal-newline mode, so a
change of line endings alone is not a change of install: numpy's, scipy's,
mpmath's and SymPy's ``RECORD`` files are written with CRLF line endings,
and the text read holds none (`[M]` 2026-10-16, the archivist's probe on
numpy 2.4.4, scipy 1.17.1, mpmath 1.3.0 and SymPy 1.14.0).

A file that several ``RECORD`` files list is refused as
:class:`~orpheus.numerics.traced_memo.Unpinnable`: pinning one owner's
version would leave an upgrade of the other unseen, a stale hit. qa's
round-two review of #602 (finding F2) showed it on a fake site directory
in which two distributions list ``ns/x.py`` and the last one read was
pinned; `[M]` 0 of the 30 446 paths in the environment are listed twice,
so no entry was wrong. Gate D2b; battery arm D-f, the refusal removed,
reddens it alone. The version is resolved by the distribution's name, not
by the ``.dist-info`` directory the table found, and the two routes agree
in the environment (`[M]` the same review, finding F3: 0 of 185
name-to-version mismatches and 0 duplicate normalised names).

Over the 474 site-packages files a
characteristic import loads, the new owner equals the old on 473; the one
that moved is ``sphinxcontrib/tikz.py`` (``scratch/characteristic_architecture/p3/m602/fix/d1_probe.py``).
The gates are ``tests/gates/numerics/test_traced_memo_ownership.py`` (D1
to D8 and D2b; D1 checks each owner against the standard library's own ``RECORD``
reader, an independent route).

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

**Every input at its first touch.** A run has three kinds of file-system
input: a file it opened for reading (by its bytes), a path it probed (by
what was there: a file, a directory, something else or nothing), and a
directory it listed (by its sorted names). The three kinds are one
table in :mod:`orpheus.numerics._traced_memo_boot`: ``Observation(field,
observe, changed_by)``, with the values ``CONTENT`` (the field
``opened``, observed by ``file_digest``), ``PRESENCE`` (``probed``,
``presence``) and ``LISTING`` (``listed``, ``listing_digest``);
``changed_by`` says whether the run's own writes change what is seen at
a path. Each input is observed at the moment the run first touches it,
by one function, ``_observe``, which keeps the observation in the
recording field the kind names (each maps a path to what was seen
there); the manifest pins those observations. Until #602's cut E all three were
observed when the manifest was written, after the run, so an input
changed while a generation ran was pinned at its new state beside an
answer computed from the old one: the entry validated and served the
stale answer (ERR-116 in :ref:`theory-verification-error-catalog`, a
time-of-check against time-of-use window; `[M]` 2026-10-16, qa's review:
1.0 served after a data file became 2.0; qa's third review, finding F2:
the same window for a presence flag and a directory listing). Observed at
the first touch, the entry is born ``Stale``, the call raises "does not
validate", and the next call regenerates from the new state. Gate
``test_qa2_a_data_file_rewritten_during_the_run_leaves_the_entry_stale``
(fresh and forked) carries ``@pytest.mark.catches("ERR-116")``. A file
rewritten and restored within one run is not seen; the run then read two
contents under one pin, which no observation of a path can name.

**An open is decoded from its flags.** The audit event ``open``
(``_OPEN_EVENTS``) carries the path, a mode and the ``O_*`` flags; both
CPython emitters (``io.open`` and ``os.open``) carry the flags, and
``io.open``'s mode string is redundant with them (``"r+"`` is
``O_RDWR``, ``"w+"`` adds ``O_TRUNC``, ``"a+"`` adds ``O_APPEND``).
``_access`` reads two answers from the flags. An open READS the path's
prior bytes when its access mode is not ``O_WRONLY`` and neither
``O_TRUNC`` nor ``O_EXCL`` is set (a truncating or exclusive open reads
nothing that was there). It WRITES when its access mode is not
``O_RDONLY``, or when ``O_CREAT``, ``O_TRUNC`` or ``O_APPEND`` is set. An
open that does both is observed first, as an input, and noted as written
second, so the refusal below applies to it when the run changed the
file. Until the elegance review of cut E an open was classified by its
mode string alone, as a read or as a write, never both, so an ``"r+"``
open (and, in E's review rounds before the fix, an ``os.open`` with
``O_RDWR``) was recorded as a write only and what it read was never
pinned (`[M]` 2026-10-16, the elegance-enforcer's
``scratch/characteristic_architecture/p3/m602/elegance4/probe_rplus.py``:
0 data pins); that hole predates E and is part of ERR-116.

**The run's own output is not an input.** The recording also holds every
path the run wrote, truncated, renamed or removed (``Recording.written``,
through ``_note_written``): an open that writes, by the rule above, and
the audit events ``os.remove``, ``os.rename`` (both paths) and
``os.truncate``. A path the
run touches only after it has written it is its own output and is not
observed: a generator that writes a scratch file, reads it back and
deletes it depends on nothing outside itself (finding N1 of qa's second
review: observed, it was born stale on every miss). And a digest reads
only a regular file (``file_digest``): a directory, a device or a pipe
digests as ``"absent"`` and is never read, which is what the file check
had always said of it (finding N3: reading ``/dev/urandom`` to its end
hung the miss).

**A run that changes an input it observed is refused.** When the run
observed a path and THEN wrote, truncated, renamed or removed it, the
state it saw no longer exists, so no state a later lookup can find pins
the answer. The manifest raises
:class:`~orpheus.numerics.traced_memo.Unpinnable`
(``_unchanged_by_the_run``, applied to every kind: a path is refused when
its kind's ``changed_by`` says the run's writes reach it and its state
now differs from the observation; the message reads "the run observed it
and then changed it"). For a file and a presence flag the writes reach
the path when the path itself is among them; for a listing, also when the
run wrote, created, renamed or removed a path directly inside the listed
directory. qa's third review, finding F1 (`[M]` 2026-10-16):
a generator that read a counter ``n`` and wrote ``n + 1`` was served as a
``Hit`` returning ``n`` under the round-three rule, which pinned a file the
run wrote by what the run left, while a cold run returned ``n + 1``. That
rule is retired. The refusal also covers a cache file read while absent
and then created (finding N2 of the second review): such a generator is
refused, never cached, a deliberate choice; no generator in the tree
uses that pattern (`[R]` the orchestrator's statement, 2026-10-16, not
re-measured for this chapter). A
directory listed and then written into is refused the same way, since
its names are then not the ones the run saw.
Gates ``test_n1_n2_a_runs_own_files_are_pinned_by_what_it_left``,
``test_n3_a_character_device_opened_for_reading_is_never_digested`` and
``test_writes_through_os_open_and_a_rename_are_the_runs_own``, on both
routes, in ``tests/gates/numerics/test_traced_memo_template.py``. The
first row was written for the retired pin-at-the-end rule and still
asserts that N2's cache file is a ``Hit``; its revision for the refusal,
and the rows for a counter, a presence flag, a listing and a truncation,
are owed by the test-architect.

**Presence probes.** While a recording is active, ``os.stat`` and
``os.lstat`` are replaced by probes that record the path, so
``os.path.exists`` and its relatives are dependencies: each probed path
becomes a :class:`~orpheus.numerics.traced_memo.PresencePin` holding what
was there at the first probe (the previous paragraph). The probes are removed when the recording stops. In P3, 0 of the 3 clients read a nuclear-data file (`[M]` the
AST census: 4 files under ``orpheus/`` read data, all in ``micro_xs`` and
the Sood cache), so the hook is for P4's generators.

**Directory listings** (``os.listdir``, ``os.scandir``) are recorded,
with the digest of the directory's sorted names at the first listing.

**The import system is exempt, decided by the immediate caller.** The
import system's own reads, probes and listings are how it finds a module,
and the module it finds is pinned where its code runs; recording them
would let any new file in any ``sys.path`` directory invalidate every entry
(`[M]` 2026-10-04: the M4.4 real-tree witnesses fail that way without the
exemption). The import system is the frozen bootstrap and the
``importlib`` package (``importlib.metadata`` lists every ``sys.path``
directory to find the installed distributions, which are pinned by
version and ``RECORD``). The exemption is decided by the frame that made the call, never
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
:func:`~orpheus.values.content.constructor_arguments` the write path
runs; an edit to any of them makes every entry stale. In a generating
process the pin is the digest of the bytes the pinned-source loader
compiled those three files from (``_RUNNING_MEMO_DIGEST``, computed once
at the memo's import from the bytes the loader recorded), never a read of
the files when the manifest is written. A read after the run is the
defect class of ERR-115: an edit to the memo during a generation would be
pinned beside an entry the old memo wrote, and the entry would validate
against the edited memo (the elegance review of #602 round two, finding
F4). Gate B15 edits the memo's source on a copy of the tree while a
generation runs and requires the entry to read ``Stale``; `[M]` battery
arm B-s, the pin read again after the run, reddens B15 alone. The same
digest names the source indexes a generation writes
(:ref:`verification-reference-cache-index-write-rule`). Validation
compares the pin with the three files on disk now, as for every other
source pin. :func:`~orpheus.numerics.traced_memo.trace_call` has no
loader and keeps the read after the call
(:ref:`verification-reference-cache-unseen`).

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
imports nothing and runs nothing (M2.4 counts it): a def pin reads its
def's digest from the source index of the file's current bytes; a module
pin reads the skeleton digest from the same index; a distribution pin
reads the installed version and the digest of its ``RECORD``; a child pin looks the child up in
the same store, recursively. Each read of a file or of what is installed is
kept for the process while every file and directory it looked at keeps its
identity on disk, so a second validation reads again only what changed
(:ref:`verification-reference-cache-kept-reads`). An index is found by the file's BYTES, so an
edited file is a new key and never a stale read (M2.5). With the indexes
on disk, a process's first validation costs 0.036 s where parsing every
pinned file cost 0.27 s, and a later one 11 to 13 ms (`[M]` 2026-10-15,
the read-back validation of the cheap-miss fixture at ``155fbdbd`` and at
``e666e58d``; :ref:`verification-reference-cache-source-index`).

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

.. _verification-reference-cache-kept-reads:

What a validation reads is kept for the process
-----------------------------------------------

**The problem, measured.** `[M]` 2026-10-17 (#602; the profile
``scratch/characteristic_architecture/p3/c3/retime/warm.prof``, read by
``an.py`` beside it). At ``main`` ``3ff7866b`` the 431 not-slow rows of
the system, reading and assembly gates ran warm in 222 s and with no memo
in 221 s: a warm store saved nothing, where on 2026-10-14 the prototype's
warm pass had saved 140 s. A cProfile of the warm run, 292.5 s under the
profiler, put 111.5 s in the 5 346 calls of
:func:`~orpheus.numerics.traced_memo.validate` (5 347 lookups). Over the
whole profiled process there were 832 287 file opens, 711 872 of them
through ``Path.read_bytes``, and 1 406 174 ``stat`` calls; the lookup of
a pinned source's file on ``sys.path`` (``_locate``) ran 690 354 times,
the first-party text reads of a validation (``_Checkout.text``) 2 282 208
times for 34.9 s, and the lookup of what is installed
(``_Installed.installed``) 38 142 times for 47.7 s.

Two facts multiplied the reads. Every validation built a fresh
``_Checkout``, and its texts and its installed environment lived for that
one validation only, so a second validation in the same process read
every pinned file and found every pinned distribution again. And an entry
pins many files: the elegance review's census over 400 entries of the
store counted, per entry, 226.2 def pins, 125.7 module pins, 62.0
presence pins, 17.5 data pins, 6.3 distribution pins, 1.0 memo pin and
0.7 child pins (`[M]` 2026-10-17). Since the transport block became a
client (``9178ad00``) a solve's blocks are its child entries, and a child
pin validates its child recursively with a checkout of its own, so one
hit of a solve read the same sources once more per child.

**The primitive.** Three objects in the standard-library-only boot module
:mod:`orpheus.numerics._traced_memo_boot` (the user's ruling of
2026-10-17: the change-detection layer is a primitive of that module, so
it moves with the harness package of #501 and imports nothing from
``orpheus``; it is the first stage of the change-detection layer of
``.claude/plans/verification_tree_refinement.md``, ruling 15):

* :class:`~orpheus.numerics._traced_memo_boot.FileIdentity` is a path's
  identity on disk as ``stat`` reads it: its device, its inode, its size,
  and its modification and change times in nanoseconds.
* :class:`~orpheus.numerics._traced_memo_boot.Witness` records the stats
  one read makes, each taken BEFORE the read of what it guards, and
  whether every one of them was settled. A stat that finds nothing is
  recorded too, as ``None``.
* :class:`~orpheus.numerics._traced_memo_boot.Kept` is a read kept for
  the life of the process. It is built from a function
  ``read(key, witness)``. Its key holds the read's inputs that are not
  files (a ``sys.path``, a working directory, a name); its witness holds
  the files and directories the read looked at. A call serves the kept
  value while a fresh ``stat`` of every witnessed path shows the
  identity recorded, a missing path still missing included, and reads
  again otherwise.

The argument that an unchanged identity means unchanged content holds on
POSIX. Any write to a file moves its change time, and no call sets the
change time back: ``os.utime`` sets the access and modification times
and moves the change time to now. A replacement by rename gives the path
another inode. A name added to, removed from or renamed in a directory
moves the directory's modification time, so a directory's identity
witnesses its names. So the identity is a cheap face of each observation
the boot module defines (a file's bytes, a directory's names), and an
unchanged identity implies an unchanged digest. Gate M2.3d makes the
change time load-bearing: an edit of the same size whose modification
time is set back with ``os.utime`` keeps the size, the modification time
and the inode, and moves only the change time.

**Stat before read.** The order is the race argument. A write that lands
after the stat and before or during the read is kept under the identity
from before the write, so the next stat sees another identity and the
read is made again. The other order, read and then stat, would file the
bytes from before a write under the identity after it, and serve them for
the life of the process. Gate M2.3g writes the file right after the read
and requires the next validation to name the edit.

**Only settled reads are kept.** A second write within one timestamp tick
can leave the size and both times as the first write left them; git calls
such index entries "racily clean". A stat is settled when the later of
its modification and change times lies more than
:data:`~orpheus.numerics._traced_memo_boot.TRUSTED_AGE_NS` (2 s, which
covers FAT's two-second timestamps, the coarsest of a common filesystem)
before the process's wall clock. A read with any unsettled stat is
returned and not kept, so the next call reads again. A read that depended
on something no stat witnesses is marked ``unkeepable`` and never kept.
Gate M2.3e simulates every write landing in one tick (an identity whose
two times never move) and requires a same-size edit of a freshly written
file to be named.

**The five kept reads**, each a module-level ``Kept`` of
:mod:`orpheus.numerics.traced_memo`:

.. list-table:: What a validation keeps for the process
   :header-rows: 1
   :widths: 16 26 36 22

   * - Kept read
     - Key
     - Witness
     - Read by
   * - ``_SOURCES``
     - the pin's path, ``tuple(sys.path)``, the working directory
     - every candidate up to the one read, a missing one included: the
       absolute path itself, or the pin's path under each ``sys.path``
       entry in order (an empty entry is the working directory)
     - ``_Checkout.text``, so every def pin and module pin
   * - ``_DISTRIBUTIONS``
     - the distribution's name, ``tuple(sys.path)``, the working
       directory
     - every ``sys.path`` directory, before the find (an install or an
       uninstall adds or removes a ``.dist-info`` there); then the found
       directory's ``METADATA``, ``PKG-INFO`` and ``RECORD``, named from
       the directory. A find that names no directory is unkeepable.
     - ``_Installed.installed``, so every distribution pin and the
       manifest's ``_distribution_pin``
   * - ``_FILE_DIGESTS``
     - the data pin's path (absolute, as the recording observed it)
     - the file
     - :class:`~orpheus.numerics.traced_memo.DataPin`
   * - ``_LISTINGS``
     - the listing pin's path
     - the directory
     - :class:`~orpheus.numerics.traced_memo.ListingPin`
   * - ``_MEMO_DIGESTS``
     - the memo's three source files (``traced_memo.py``,
       ``_traced_memo_boot.py``, ``content.py``)
     - the three files
     - :class:`~orpheus.numerics.traced_memo.MemoPin`, and the memo pin
       of a manifest built outside a generation (``trace_call``'s)

``sys.path`` and the working directory are in the KEY, never in the
witness, because they are not files: a different ``sys.path`` is a
different read, made fresh. Shadowing is why ``_SOURCES`` witnesses every
candidate and not only the file it read. A pin's path is relative to its
``sys.path`` entry (``_pin_path``), so a module can appear on an EARLIER
entry that is already on ``sys.path``; the same import would then load
it. The earlier candidate was witnessed as missing, so its appearance is
a change, the read is made again, and the shadowing file is the one read
(gate M2.3f). ``_DISTRIBUTIONS`` names the files it guards from the found
directory and not from what ``RECORD`` lists, so a distribution whose
``RECORD`` does not list its ``METADATA`` is still guarded (the elegance
review, A6). :func:`~orpheus.numerics._traced_memo_boot.file_digest` and
:func:`~orpheus.numerics._traced_memo_boot.listing_digest` take an
optional witness, so a recording's observation of a file and a pin's
validation of it are one definition.

What stays per validation: ``_Checkout._texts`` gives one state per
source file within one validation, so every pin of a file is judged
against one read of it; it does not reach across files, nor across
children, each of which validates with a checkout of its own. The
ownership table of ``_Installed`` (which distribution's ``RECORD`` lists
a file) is still built once per checkout; what is installed as each
distribution is the process's ``_DISTRIBUTIONS``. Presence pins are not
kept: a presence probe is itself one ``stat``. ``_locate`` and
``_Installed._versions`` retired.

**The scope boundary.** The argument above holds on a POSIX filesystem
whose timestamps come from the clock this process reads. Outside it:

* On Windows ``st_ctime`` is the creation time, so the same-size edit
  with its modification time set back (M2.3d's) would be served stale
  (qa's review, finding 7a).
* On a network filesystem the timestamps come from the server's clock,
  which may differ from the client's, and attribute caching may serve an
  old identity (finding 7c); a wall clock stepped backwards on a
  coarse-tick filesystem, with a same-size edit in place, is the third
  case (finding 7b).
* The racy reference is the process's wall clock. git compares against
  the index file's own modification time, so both sides of its
  comparison come from the filesystem's clock. Inside one process the
  difference is small; a persisted identity store, the next stage of the
  change-detection layer, takes its own write time as the reference (the
  elegance review, A5). `[R]`
* Beneath ``_DISTRIBUTIONS``, the standard library's
  ``importlib.metadata`` keeps each ``sys.path`` directory's listing in a
  cache keyed on the directory's float ``st_mtime``, with no racy guard
  of its own. When ``_DISTRIBUTIONS`` declines to keep a find whose
  directory was written within the trusted age, the find it repeats can
  still be handed a stale listing by that cache. The exposure predates
  this change (qa's review, finding 8; `[R]`, unmeasured).

None of these is the host on which the memo is measured (APFS, a local
disk). qa's review of the first version found no realizable state on
this host in which a kept value is served while the file on disk differs
(`[M]` 2026-10-17,
``scratch/characteristic_architecture/p3/c3/hitcost/qa/review.md``, ten
attacks withdrawn with their reasons).

**What was tried first, and what the reviews changed.** The first version
of the fix (``53551868``, never merged) had two hand-written caches,
``_Sources`` and ``_Distributions``. Each spelled the four steps itself
(stat before read; keep only if settled; serve while the stats are equal;
read again otherwise), and they spelled settledness differently: one
called a predicate inline, the other stored an ``"unsettled"`` sentinel as
a third value beside an identity and ``None``, a state the kept record's
type admitted and the code had to keep out. Data pins, listing pins and
the memo pin still read afresh, so the subject line "a validation
re-reads only what changed on disk" held for 2 of the 5 pin kinds that
read a file. The distribution's guard set was filtered from what its
``RECORD`` lists. And the ladder marker of the old M2.3b had stayed on
the new helper ``_install``, which pytest never collects. The elegance
review (A1 to A9, ``scratch/characteristic_architecture/p3/c3/hitcost/elegance_review.md``)
proposed one stage-1 primitive in the boot module: ``FileIdentity``,
``Witness``, ``Kept``, with the inputs that are not files in the key. The
landed version is that primitive, with all five kept reads built on it.

qa's review of the first version (``.../hitcost/qa/review.md``) applied
five arms, each removing one claim the commit made, to the 23 shipped
rows. Three of the five reddened no shipped row: the claims of shadowing,
of stat before read and of ``METADATA`` in the guard set. The racy
guard's row, M2.3e, reddened under its arm but skipped by chance in 2 of
12 isolated runs, because it waited on the wall clock, so a merge run
could read it as green. Those four claims became gates M2.3f, M2.3g and
M2.3h and the deterministic M2.3e. M2.3i gates the
data and listing pins the elegance review absorbed. M2.3b was re-posed:
it had faked an upgrade by monkeypatching
``importlib.metadata.distribution``, a change no disk shows, and it now
installs a numpy ``9.9.9`` into a directory already on ``sys.path`` and
then rewrites that install's ``RECORD`` in place.

**The gates and their first reds.** The rows are in
``tests/gates/numerics/test_traced_memo_validation.py``. Each arm is
applied at ``pytest_configure`` in the process pytest runs, over the 27
rows of the file (`[M]` 2026-10-17,
``scratch/characteristic_architecture/p3/c3/hitcost/arms2.log``):

.. list-table:: The kept-read gates
   :header-rows: 1
   :widths: 9 55 36

   * - Row
     - Assertion
     - First red (arm: rows red)
   * - M2.3b
     - An install into a ``sys.path`` directory that was already there,
       then its ``RECORD`` rewritten in place, are each a reason in the
       next validation of the same process.
     - ``distributions-blind``, a kept find served without its stats:
       M2.3b, M2.3h
   * - M2.3c
     - A second validation with nothing changed makes no read of a
       pinned source (the spy counts the reads ``_SOURCES`` makes,
       however a read is spelled).
     - ``never-serve``: M2.3c
   * - M2.3d
     - A same-size edit with its modification time set back is a reason.
     - ``drop-ctime``: M2.3d, M2.3i; ``always-serve``: M2.3b, d, f, g,
       h, i
   * - M2.3e
     - A write within the tick of the last read is read again
       (deterministic: an identity whose two times never move).
     - ``settled-true``, the racy guard removed: M2.3e
   * - M2.3f
     - A module appearing on an earlier ``sys.path`` entry that shadows
       the pinned file is a reason.
     - ``witness-the-read-file-only``: M2.3f
   * - M2.3g
     - A write between the read and the next stat is read again.
     - ``read-then-stat``: M2.3g
   * - M2.3h
     - A ``METADATA`` rewritten in place, same size, is a reason.
     - ``record-only``: M2.3h
   * - M2.3i
     - A data pin's file edited in place with its time set back, and a
       name added to a listing pin's directory, are reasons.
     - ``never-witness``, a kept read served without its stats: M2.3i

The positive control, ``control-empty-sources`` (every source read
returns nothing), reddens 21 of the 27 rows.

**The measured outcome.** The three files' timing, #602's last done-when
("the three-file cold/warm/no-memo timing is re-measured and recorded on
the reference-cache theory page"), is the table below. Each row is one
serial run under ``-O`` with ``OMP_NUM_THREADS=1`` in a fresh worktree
(an empty store) of the tree named, over the same 431 not-slow rows
(`[M]`; the drivers are ``run.sh`` beside each summary):

.. list-table:: The system, reading and assembly gates (not slow), seconds of wall time
   :header-rows: 1
   :widths: 30 12 12 12 34

   * - Tree
     - No memo
     - Cold
     - Warm
     - Source
   * - The block memo's prototype, 2026-10-14
     - 242
     - 1 303
     - 102
     - ``p3/c1/gate_spec.md`` §7; 383 entries
   * - ``155fbdbd``, #602 round two, 2026-10-15
     - 236
     - 406
     - 168
     - "What a miss costs" below; 257 block entries
   * - ``main`` ``3ff7866b``, 2026-10-17
     - 221
     - 816
     - 222
     - ``p3/c3/retime/summary.txt``; 1 123 entries
   * - ``53551868``, the first version of this fix, 2026-10-17
     - 227
     - 725
     - 155
     - ``p3/c3/hitcost/summary.txt``; 1 123 entries
   * - ``8fbb1d61``, the landed fix, 2026-10-17
     - 221
     - 820
     - 163
     - ``p3/c3/hitcost2/summary.txt``; 1 123 entries

The paths are under ``scratch/characteristic_architecture/``. Every row
is a single run. Other agents' work shared the host during the cold and
warm runs of the ``53551868`` row; the last row ran alone. Each
``summary.txt`` records the ``HEAD`` it timed. The last row was timed at
``b57648b5``, whose tree is ``8fbb1d61``'s (the two differ in the commit
message only). The landed version adds to the first the shadowing
witness, the data, listing and memo kept reads and the directory-named
guard set. Its warm pass, 163 s, is the first version's 155 s within what
one run of each can separate; its cold pass, 820 s, is ``3ff7866b``'s
816 s, so the first version's 725 s cold is `[R]` host variation rather
than a saving of the fix: a miss re-reads nothing the fix keeps, since a
miss validates nothing.

Between ``155fbdbd`` and ``3ff7866b`` the store grew from 257 block
entries to 1 123 entries of every kind, as more of the readings' rows
read through the memo (``6af57441``, ``410bfca8``); `[R]` that, with
the per-validation re-reads, is why the warm pass rose from 168 s to
222 s. A warm pass of the fix costs 67 s less than ``3ff7866b``'s and
72 s less than no memo. It is still above the prototype's 102 s, for
two reasons `[R]`, neither profiled after the fix. Each entry now pins
about 130 modules, where the prototype's blocks pinned 13
(``gate_spec.md`` §5), and the store holds 1 123 entries against 383;
a kept read still makes one ``stat`` per witnessed path, which for a
relative source pin is one per ``sys.path`` entry up to the one that
holds the file, and every pin is still walked. And part of the run is
computation outside the memo: the assembly gates build their blocks in
process, and the reading gates' point rows and ``bypass()`` rows are
computed (the paragraph after the round-two table below).


.. _verification-reference-cache-source-index:

The source index: each source text parsed once per content
==========================================================

A manifest and a validation read the same few facts from a first-party
source file: the digest of its skeleton (the module pin), the digest of
each outermost def that ran (the def pins), and, to map a code object
that ran onto the def that holds it, the span of every node a code
object can be compiled from. Those facts are a function of two things
only: the file's bytes, and the memo code that normalises and digests
them. A fresh process holds none of them in memory, so without a store
every miss parsed each of the about 128 first-party files its imports
touched, and a process's first validation parsed every file its entries
pinned. The source index keeps the facts on disk under the content that
determines them: the first generation that meets a text parses it once
and writes its index, and every later process reads the index back.

What an index holds, and why it holds no tree
---------------------------------------------

An index is a ``_SourceIndex`` (a frozen value, private to
:mod:`orpheus.numerics.traced_memo`), built from one parse of a text's
bytes by ``_parsed_index(data)``. It has three parts:

* ``skeleton``: the digest of the module's skeleton, the syntax tree with
  every outermost def's body read as ``pass`` and every docstring removed
  (:ref:`verification-reference-cache-manifest`).
* ``defs``: one ``_IndexedDef`` per outermost def, in the order of the
  walk that finds them. Each holds the def's qualified name (``Class.method``
  for a def in a class body), its ordinal among the defs of that name (0
  for the first), the first and last line of its span (inclusive; the
  first line is that of its first decorator), ``body_line`` (the line of
  the first statement of its body once the docstring is stripped) and the
  digest of its normalised source.
* ``scopes``: one ``_Scope`` per node that a code object is compiled
  from, anywhere in the file (a ``def``, an ``async def``, a ``class``, a
  ``lambda``, a generator expression, and a list, set or dict
  comprehension). Each holds the name of that code object (the def's or
  the class's own name, or ``<lambda>``, ``<genexpr>``, ``<listcomp>``,
  ``<setcomp>``, ``<dictcomp>``) and the node's inclusive span.

Three queries read it:

* ``digest(qualname, ordinal)``: the digest of the ``ordinal``-th def of
  that name (``-1``, the default, is the last, the one Python binds), or
  ``None`` when there is none, which ``DefPin.stale`` reports as "gone".
  The public :func:`~orpheus.numerics.traced_memo.function_digest` is this
  query, and :func:`~orpheus.numerics.traced_memo.skeleton_digest` reads
  ``skeleton``.
* ``enclosing(line)``: the outermost def whose span holds ``line``, or
  ``None``.
* ``holds(within, line, name)``: whether a node whose code object is named
  ``name`` spans ``line`` and lies inside the span of the def ``within``
  (``None``: anywhere in the module). The name of a signature scope
  (``__annotate__``, a type-parameter scope) answers yes: the skeleton
  pins its source.

``_def_pin`` composes the three. A code object first reported at ``line``
under the name ``name`` is module-level code, which the skeleton pins
(no def pin), when it is ``<module>`` or a signature scope that starts
before its def's ``body_line``. Otherwise the def that holds it is
``enclosing(line)``, and ``holds`` must confirm that a node of that name
spans the line inside that def; the pin is that def's ``DefPin``, or none
when no def encloses the line. A code object no node accounts for is
refused as :class:`~orpheus.numerics.traced_memo.Unpinnable`
(:ref:`verification-reference-cache-manifest`, "Code that is not a def
body").

**It holds no tree**, for three reasons.

* It is written to disk and read back, and the facts above are a few
  kilobytes of JSON per file where a syntax tree has no serial form short
  of the source it was parsed from.
* The in-memory cache keeps only the index. ``_parsed_index`` is cached
  by the text's bytes (an edited file is another key, never a stale read),
  and caching the index rather than the parsed tree does not keep every
  module's tree alive for the life of the process.
* One source text has one representation. Before ``155fbdbd`` the parse
  was a second type beside the index, a parsed text holding its tree and
  its defs, and the two answered ``digest`` with α-identical bodies (the
  elegance review of #602 round two, finding F1): one concept twice, which
  would let the public :func:`~orpheus.numerics.traced_memo.function_digest`
  and the pins drift apart under a change to one of them. The parse is now
  the index's constructor, and one type answers every query.

**Every def is digested when the index is built.** A miss that parses
each file pins only the defs that ran (247 of 1 948 on the cheap-miss
fixture), so digesting a def only when a pin asks for it saved work while
every miss parsed afresh. An index written to disk cannot know which defs
a later process's pins will ask for, and it is built once per content and
read many times, so it digests them all. `[M]` 2026-10-16, this chapter's
probe (``scratch/characteristic_architecture/p3/m602/archivist_numbers.py``) over the 390 tracked ``orpheus/*.py`` files, ``-O``, in memory:
building their indexes from their bytes takes 1.40 s (4 198 defs, 6 366
scopes, 713 kB of JSON, the largest index 35 kB), and reading the 390
back from disk, with every header and body digest checked, takes 33 ms
(median of five passes). Over the 128 files one characteristic block
imports, the specification measured a build of 557 ms, a read-back of
9.3 to 9.6 ms (minimum and median of seven) and a cold parse of 225 ms
(``scratch/characteristic_architecture/p3/m602/round2_gate_spec.md``
§1.1, on the draft format).

**Why** ``holds`` **asks span containment.** Outermost defs never
overlap: the walk that finds them never descends a def's body, and a def
is a statement, never part of an expression. So a node lies inside an
outermost def's span exactly when it lies in that def's syntax tree, and
the question "is this node part of the def that ran" is two comparisons
of line numbers, which the index can answer without the tree. The
specification's draft stored, beside every scope, the position of the def
that owns it, with ``-1`` for none; the elegance review found the field
derivable and the sentinel dangerous. `[M]` 2026-10-15, its probes
``scratch/characteristic_architecture/p3/m602/elegance2/probe_owner.py``
and ``probe_control.py`` over the tracked files: the stored owner equals
the def whose span contains the scope's span on 24 262 of 24 262 scopes,
and shifting the owner of each of the 5 664 scopes that have one makes
5 664 of 5 664 mismatch (the control). A ``-1`` read as a list index
names the last def, silently, and ``-1`` meant "no owner" in one place and
"any owner" in another. The field and the sentinel are gone; ``enclosing``
returns the def or ``None``.

**The key** names an index by the digest of three components, joined by
a NUL byte:

* the digest of the text's **bytes**: an edited file is another name, so
  an index never goes stale, it is only no longer asked for (gate B6);
* the **interpreter** (``sys.version`` and the bytecode cache tag): the
  ``ast.dump`` that the digests hash is not a stable interface across
  Python versions, although `[M]` 0 of 390 ``orpheus/`` files differ in
  it between 3.13.12 and 3.14.3 (the specification's
  ``ta/p2_ast_across_pythons.py``), so the component is mechanical (gate
  B7);
* the **memo's source**: the docstring strip, the walk and the dump are
  code in the memo, a change to them changes every digest, and naming the
  index by the memo's source makes that change a change of name (gate
  B8). Which digest of the memo's source this is, is the write rule below.

The ``-O`` level is not part of the key: the build contains no ``assert``
and ``ast.parse`` does not optimise, and gate B9 shows indexes built with
and without ``-O`` are byte-equal.

**The file.** An index lives at ``.cache/sources/<name>.json`` under the
repository root, beside ``references`` (the entries) and ``bytecode`` (the
pinned-source loader's cache). The three stores share one root, ``CACHE``,
defined once in :mod:`orpheus.numerics._traced_memo_boot` as the real path
of the repository's ``.cache``; before ``155fbdbd`` the root was spelled
three times, two of them resolving symbolic links and one not (the
elegance review, finding F3). The file is one header line, then the body:

.. code-block:: text

   {"bytes": <digest of the text>, "python": <interpreter>, "memo": <digest of the memo's source>, "body": <digest of the body>}
   {"skeleton": <digest>, "defs": [[qualname, ordinal, first, last, body_line, digest], ...], "scopes": [[name, first, last], ...]}

**Reading.** ``_source_index(data)`` serves an index from memory (a
cache by the bytes), else from its file, else from a parse. ``_read_index``
returns the index, nothing (no file), or a ``_CorruptIndex`` with its
reason. An index is corrupt when its header names other key parts (the
index of another text copied under this name), when its body does not
match the body digest in its header (a flipped bit), when it does not
parse into the three parts (a truncated or malformed file), or when the
file cannot be read at all (an ``OSError`` other than a missing file, a
``PermissionError`` for one; before ``155fbdbd`` that error escaped and
crashed the lookup, qa's round-two finding F4). A corrupt index is never
read: the text is parsed and its index served, and a generating process
writes the file again over it. The reads run under
``_traced_memo_boot.unrecorded``: an index is the memo's own reading, so
no recording pins a file under ``.cache/sources`` (gates B14 and B18).

**Writing** goes through ``write_atomically`` in
:mod:`orpheus.numerics._traced_memo_boot`, the one atomic write of the
memo's three stores: the bytes go to ``<path>.<pid>.<thread id>.part`` in
the same directory and are moved into place with one ``os.replace``, so a
reader finds a whole index or none (gate B11, four processes writing and
reading one name 200 times each). Before ``155fbdbd`` the three stores
had three spellings of that write with three temporary-name schemes (the
elegance review, finding F3).

.. _verification-reference-cache-index-write-rule:

Only a generating process writes an index
-----------------------------------------

The name of an index is meant to name the code that built it. A
generating process can name that code exactly; no other process can.

**A generating process knows the bytes it ran.** The pinned-source loader
(:ref:`verification-reference-cache-pinned-source-loader`) compiled the
memo's three source files, ``traced_memo.py``, ``_traced_memo_boot.py``
and ``content.py``, from bytes it read and recorded in
``_traced_memo_boot.PINNED_SOURCES``. At the memo's import,
``_RUNNING_MEMO_DIGEST`` is the digest of those recorded bytes: the code
that is running. A generating process names every index it writes by it,
and the entry it writes pins the memo by it too (``MemoPin``, the
manifest section).

**Any other process cannot.** The test process, or a script, imported the
memo through the ordinary loader, which may have run a timestamp ``.pyc``
of a text since edited at the same size within the same second (the
mechanism of ERR-115, applied to the memo itself), and the source on disk
may change between the loader's read and any later read of it. Such a
process has ``_RUNNING_MEMO_DIGEST = None``. It names indexes by
``_INDEX_MEMO_DIGEST``, the digest of the memo's source as it stood on
disk when the memo was imported; it reads the indexes that generations
wrote under that name, parses in memory when none exists, and never
writes. Its read of the memo's source at import is nobody's dependency,
since no recording is active when a process that is not generating
imports the memo.

**The founding case** (qa's round-two review of #602, finding F1; `[M]`
2026-10-15, ``scratch/characteristic_architecture/p3/m602/qa2/p2_stale_pyc_poison.py``
and its log, on a temporary copy of the memo's files). Before the rule,
any process wrote an index, named by the memo's source as read from disk.
Process A imports the memo, whose normaliser strips docstrings. A
same-size edit makes the normaliser keep them, and the file's
modification time is put back. Process B imports the OLD normaliser from
the ``.pyc``, names its index by the NEW source, and writes the index of
a sample file. Process C compiles the new code afresh and reads B's index
under its own name: the index gives the digest ``c923c2ad`` where C's own
parse gives ``e3c45572``. The control, C with the index deleted, gives
``e3c45572`` both ways. No check on the file could see it, since a
well-formed index is accepted whoever wrote it, and an index is never
re-validated, so the wrong digests stand for that text as long as the
file exists. A stale hit then needs a coarser poisoned normaliser that
gives the pinned text and an edited text one digest `[R]`; the poisoning
is `[M]`.

**Why the indexes a non-generating process reads are the right ones**
`[R]`. Such a process reads indexes to validate entries. An entry is a
``Hit`` only if its ``MemoPin`` equals the memo's source on disk now, and
that pin is the code of the generation that wrote the entry and its
indexes. So an entry that can be served was pinned by digests built by
the code whose name the reading process looks indexes up under, and the
validation compares like with like. An index written by a memo that is
not the one on disk is never read by a validation that could serve.

The gates: B16 (both sides of the rule: outside a generation an index is
parsed in memory and no file is written, and in the generating condition
the same call writes one) and B17 (qa's scenario on a copy of the memo's
files: C is served its own normalisation's digest). `[M]` battery arm B-q,
which lets every process write, reddens B16 and B17.

``trace_call``, the in-process instrument, runs in a process that is not
generating, so it writes no index, and its memo pin is the digest of the
memo's source on disk after the call, read through the process's kept
read ``_MEMO_DIGESTS``, which reads the three files again when any of them
has changed (:ref:`verification-reference-cache-unseen`).

.. rubric:: Gotcha: a generation's manifest depends on the state of the index directory

A generation that validates a child entry runs the memo's index code
while it records, and ``unrecorded`` silences the audit events and the
probes but not ``sys.monitoring``, which reports every code object that
starts. A cold index directory runs the build (``_parsed_index``, the
JSON writer) and a warm one runs the reader, so the def pins of
``traced_memo.py``'s own defs differ between a cold and a warm run of the
same generation (qa's round-two review, an observation `[R]`). The extra
pins only ever make an entry stricter, never a stale hit; gate B5, which
requires one manifest cold and warm, uses five blocks with no nested
lookup.

The gates and their witnesses
-----------------------------

The rows are ``tests/gates/numerics/test_traced_memo_source_index.py``.
B1 to B4 compare the index production builds with a parse the test makes
for itself (the docstring strip, the outermost-def walk, the skeleton and
the digests re-spelled from their definitions, and, for B3 and B4, the
retired tree-walking ``_def_pin`` and ``_holds`` of ``main`` before
``155fbdbd`` as the oracle), so production's builder is on neither side;
only ``ast`` and ``hashlib`` are shared. The battery
(``scratch/characteristic_architecture/p3/m602/ta/batteryR/``, the
test-architect, 2026-10-15) applies each arm as a textual edit to a copy
of the tree and runs the 174 traced-memo rows in that copy: 174 of 174
pass unmutated, and the positive control, arm B-c (the key without the
bytes), reddens 113 of them.

.. list-table:: The source-index gates
   :header-rows: 1
   :widths: 9 55 36

   * - Row
     - Assertion
     - First red (battery arm, rows red)
   * - B1
     - An index survives its JSON round trip, over every tracked ``.py``
       file and a fixture of 26 def shapes.
     - ``from_json`` and ``to_json`` disagreeing (not in the battery)
   * - B2
     - The skeleton digest and every def digest equal the test's own
       parse.
     - A normalisation that differs from its definition (not in the
       battery)
   * - B3
     - ``_def_pin`` read from the index equals the tree-walking oracle on
       every code object reachable from ``compile`` of the file.
     - B-k, spans written exclusive (B3, B4)
   * - B4
     - The same at the first and last line of every span and one line
       outside it, which B3's realisable code objects never probe.
     - B-d, the containment check removed (B4); B-k
   * - B5
     - A real manifest of five blocks is the same from an empty and from
       a warm index directory.
     - B-m, every index read made a miss (13 rows, B5, B6 and B10 among
       them)
   * - B6, B7, B8
     - An edit to the bytes, another interpreter identity, another
       normaliser: each is another name, and the old index is not read.
     - B-c (113 rows), B-a, B-b (one row each)
   * - B9
     - Indexes built with and without ``-O`` are byte-equal.
     - None expected: a ``-O`` dependence brought into the build
   * - B10
     - A truncated (0, 1, half, all but one byte), bit-flipped, copied or
       re-headed index is corrupt, the parse's digests are served, and
       the file is valid again afterwards.
     - B-f, no body digest check (the flipped digest); B-g, no header
       check (the copied index, the re-headed one)
   * - B11
     - Concurrent writers never expose a partial index.
     - A write in place (not in the round-two battery)
   * - B12, B13
     - A process's first validation reads indexes and parses nothing; a
       miss indexes only bytes no generation has indexed.
     - The index not consulted (not in the battery)
   * - B14, B18
     - Neither the index directory nor the memo's own source is a
       dependency of a run, including a recorded run that validates a
       manifest itself.
     - B-o, the read outside ``unrecorded`` (B18)
   * - B15
     - The memo's own source edited during a generation leaves the entry
       ``Stale``.
     - B-s, the memo pin re-read after the run
   * - B16, B17
     - The write rule above.
     - B-q, any process writes (both)
   * - B19
     - An unreadable index is corrupt, and the parse is served.
     - B-r, an ``OSError`` escaping the read

One arm was blind and is removed: B-n, a quiet scope around the import-
time read of the memo's digest, which no recording can be active for, so
no row could see it. Arm B-i (an index cached by the file's path, not its
bytes) reddens 20 rows across the memo's gates, M2.5 among them.

**Equivalence with the tree route** `[M]` 2026-10-15
(``scratch/characteristic_architecture/p3/m602/fix/b3_probe.py``, the
commit's probe): ``_def_pin`` read from the index equals ``main``'s
tree-based ``_def_pin`` on 33 750 code objects of 1 236 files, 0
mismatches, with the indexes built cold and read back from disk. qa's
review: 128 of 128 indexes on disk equal a fresh build (control: a
perturbed index compares unequal), and 10 of 10 corrupt shapes read as
corrupt.


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
store's root, the chain of calls the generating processes are already
answering, and the environment the generation is handed), pickled through
the constructor pickler and sent, as one length-prefixed frame after an
envelope frame carrying the parent's ``sys.path``, on the standard input of

.. code-block:: text

   <the parent's interpreter> [-O…] -P orpheus/numerics/_traced_memo_boot.py

started in the parent's working directory with a declared environment only
(below). That is the fresh route; on a platform where forking is measured
the same command with ``--template`` starts a template, which reads one
job frame per miss and forks a child for each
(:ref:`verification-reference-cache-fork-template`). One value,
``_Launch``, holds the interpreter, the ``-O`` level, ``sys.path``, the
working directory and the environment, and starts either.

**The memo crosses by its name, and is found through the decorators
above it.** The job carries the memo as a reference, never as code:
``TracedMemo.__reduce__`` pickles it as a call of ``_resolve_memo`` on its
``function_id``, so the generating process imports the memo's module
(compiled from the bytes it pins, below) and finds the memo there.
``_resolve_memo`` reads the module attribute the ``function_id`` names and
follows its ``__wrapped__`` chain with :func:`inspect.unwrap`, stopping at
the first ``TracedMemo`` (``inspect.unwrap`` raises ``ValueError`` on a
cyclic chain). It refuses by name, with ``TypeError`` ("names no memo,
directly or through __wrapped__"), unless what it reached is a memo whose
``function_id`` is the one asked for. The unwrap is needed because a
decorator stacked above the memo (``functools.lru_cache``, or any wrapper
built with ``functools.wraps``) makes the module attribute the outer
wrapper, whose ``__wrapped__`` is the memo. Before ``7ca91c43``
(2026-10-17) the child resolved the attribute alone: it called the
wrapper, the wrapper called the memo with the key the child was
answering, and the cycle guard below refused it with "calls itself with
the same arguments". `[M]` 2026-10-17: the first stacked client,
``phi_cylinder_vacuum`` under ``lru_cache``, failed
``test_a_small_cylinders_layer_is_graded_outside_slow[interface]`` with
exactly that error. Two gates in
``tests/gates/numerics/test_traced_memo_process.py``, each with its first
red:

* **M3.1b**: a memo under ``lru_cache`` in the synthetic ``beta`` module
  generates once and returns its body's value. The old ``__reduce__``
  restored in process reddens M3.1b and leaves M3.1 green; qa reddened it
  on the fresh-interpreter route too, with forking off.
* **M3.1c**: a name whose attribute is no memo, or wraps a memo of another
  ``function_id``, is refused by name rather than run. Dropping the
  identity check reddens the another-memo row only.

No client stacks a decorator above its memo today: the first that did,
the oracle's ``lru_cache``, was removed in review (the clients, below,
say why). The resolution is the door for any later wrapping decorator.

**The boot script imports only the standard library**, so it can start
recording before any ``orpheus`` code runs: the shared start of a
generating process and a template, ``_boot``, starts a :class:`~orpheus.numerics._traced_memo_boot.Recording`,
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
would be offered to nobody. The boot's ``_boot`` therefore registers itself
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
``TMPDIR`` and the locale variables, each when the parent has it) and
``_FIXED_ENVIRONMENT``, the variables it receives at one value whatever the
caller set: ``PYTHONHASHSEED=0``, so that iteration over a set of strings
has one order in every generation, and the five BLAS and OpenMP thread
counts (``OMP_NUM_THREADS``, ``OPENBLAS_NUM_THREADS``,
``MKL_NUM_THREADS``, ``NUMEXPR_NUM_THREADS`` and
``VECLIB_MAXIMUM_THREADS``) at ``1``, so that a reduction is summed in one
order. The thread counts were passed through from the caller and pinned
until 2026-10-26; the store then held one population per shell (`[M]`
5593 entries written under a plain shell, 4451 under
``*_NUM_THREADS=1``), and a run in the other shell regenerated every
entry. A thread count changes the order in which BLAS sums a product, so
it belongs to the answer; fixing it, rather than dropping it from the
pin, keeps it there and makes it one value. The
pin, an :class:`~orpheus.numerics.traced_memo.EnvironmentPin`, holds the
declared variables as the run RECEIVED them, ``PATH`` excluded: ``PATH``
is passed (``uname`` is found through it) but the one program a
generation may start is admitted by its system path and pinned by its
bytes, and a pinned ``PATH`` made every entry stale when a virtual
environment was activated; and a pin taken after the run would never
validate for a generator that sets a declared variable (third review,
``test_e_…``). "As the run received them" means as the asking process
HANDED them (``_Job.environment``), not the generating interpreter's own
``os.environ`` after its start, which PEP 538 locale coercion rewrites
under ``LANG=C`` (:ref:`verification-reference-cache-fork-template`). Two reasons: no
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
nothing. The census found two such sites: the solve spy of the
``test_trajectory_resolvent_reference`` module (then under
``tests/gates/derivations/``, deleted with its family on 2026-10-10),
which patched both multi-region solvers to count lazy solves (its rows
read under an ``in_process`` fixture that entered the bypass), and a
counting spy on ``Symbolic.steps``. The characteristic reference's rows
that read ``answer`` directly enter the bypass the same way, because
``answer`` reads the memoised solve. M3.10 gates both halves: a patch in the parent
never reaches an entry, and a call under the bypass reads and writes
nothing.

**What this does to a mutation battery.** For the same reason, an
in-process mutation (a monkeypatch, or a ``-p`` plugin rebinding a
symbol) is invisible to every row that reads a value through the memo: a
miss generates from the source on disk, and a hit serves what that source
computed. Since the transport blocks became a client (2026-10-15), that
covers every characteristic row that reads ``GalerkinSystem.groups``, and
since 2026-10-17 every row that reads one of the four stored oracle routes
and the eight slow cylinder point-reading rows that read through
``CharacteristicDerivation.evaluate`` (concept 2; the clients, below). A
battery or a ``catches`` re-drop on such rows mutates textually (a copied
tree or a worktree, so the defect is on disk) or has its plugin enter
:func:`~orpheus.numerics.traced_memo.bypass` at ``pytest_configure``;
otherwise it reads 0 reds and reports a catcher as blind, or the tree as
clean (``vv-principles`` #17).

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
arise, and ``ExactInfiniteMediumDerivation`` kept its one field,
``medium``. Since 2026-10-08 it carries a second, ``gauge``, the
production per unit flux of each group that the question declares
(:ref:`structured-geometry-question-values-gauge`); its content is still
derived from the specification, and the rule above does not change.

The rule this ruling sets for later clients: **a reading cheaper than one
interpreter start plus one validation is not memoised.**


.. _verification-reference-cache-clients:

The clients
===========

The clients are the characteristic reading, its solve and its
transport blocks, four mpmath routes of a test oracle, and the PS-1982
vacuum-sphere solver, below. The
two trajectory-resolvent clients described first were P3's, and were
deleted with their family; their description is the record of how ruling
2 was first realised.

**The two multi-region solvers (deleted)**,
``solve_greens_function_sphere_mr`` (in
``trajectory_resolvent.greens_function``) and
``solve_greens_function_cylinder_mr`` (in ``greens_function_cylinder``),
were decorated ``@traced_memo(canonical=SOLVER_ARRAY_ARGUMENTS)``. The
census showed why the memo sat at the solver functions themselves: they
returned frozen
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

**The trajectory-resolvent reading (deleted).**
``TrajectoryResolventDerivation`` had content identity (its six init
fields; its ``Billiard`` and its rays class were derived from them and
excluded), and its ``evaluate`` was a memoised METHOD. A reading was
generated once, in a fresh process that constructed the derivation from
its content (finding F1), solved through the solver's memo (a child
entry) and evaluated the observable. The memo was on ``evaluate`` itself,
not on a module-level reading function beside it, which would be a
second spelling of ``evaluate`` free to drift from it; the characteristic
reading keeps that choice. Within one process the derivation held its
solve in a ``cached_property``; across processes, the reading and the
solve were memo entries (:ref:`characteristic-origins-history-reading`).

**The characteristic reading and its solve.**
:class:`~orpheus.derivations.continuous.characteristic.reference.CharacteristicDerivation`
has content identity (its two fields, the specification and the
resolution), and two of its methods are traced memos: ``evaluate``,
keyed on the derivation and the observable, and ``solve``, keyed on the
derivation alone. A reading's generation reads ``answer``, which calls
``solve``, so the solve is a child entry of every reading and all the
observables of one derivation share it: ruling 2 again, realised on a
method pair rather than on a module-level solver. Within one process the
``answer`` property holds the solve's decoded value. `[M]` 2026-10-10,
the elegance review's probe: three observables of one cold fundamental
wrote 3 reading entries and 1 solve entry. The page of record is
:ref:`characteristic-door`.

**The characteristic transport block.**
:obj:`~orpheus.derivations.continuous.characteristic.assembly.group_transport`
is a module-level traced memo, one group's block, called once per group
by ``GalerkinSystem.groups``, so the solve's generation reads each block
as a child entry and a point reading pins the same block keys. It is the
first client whose arguments are geometry: it needed
``PanelBasis``, ``ConcentricPartition`` and ``RigidMotion`` to become
:class:`~orpheus.values.content.ContentIdentity` values, and it declares
two canonical forms, ``sigma_t`` and ``support``, each the function the
body itself calls first, so a served block and one computed under
``bypass()`` (which passes the raw arguments) are one value. Its point is
the decoupling: the block's manifest pins the transport kernel and no
consumer, so an edit to the system, the reading or the door leaves the
expensive blocks valid while the cheap solve regenerates. The design, the
key's completeness, why its ``support`` is required, why the line rule's
chunk and budget are pinned through its source rather than passed, and
the measurements are on :ref:`characteristic-block-memo`. `[M]`
2026-10-14 on 431 not-slow rows of the three files that read blocks:
242 s without the memo, 1 303 s cold, 102 s warm; the cold overhead is
the miss cost (#602), and the user ruled on 2026-10-15 that the block is
memoised whatever its size, with no threshold. After #602 the cold pass
is 406 s against 236 s without the memo and 168 s warm (`[M]` 2026-10-15
at ``155fbdbd``; the section "What a miss costs"). On 2026-10-17, with
1 123 entries in the store, ``main`` at ``3ff7866b`` ran 816 s cold,
222 s warm and 221 s without the memo, and the first version of the
hit-cost fix (``53551868``) 725 s, 155 s and 227 s, and the landed fix
(``8fbb1d61``) 820 s, 163 s and 221 s (`[M]`, single runs;
:ref:`verification-reference-cache-kept-reads`).

**The mpmath routes of a test oracle** (#405 P3 concept 1, 2026-10-17,
``333268fb``). ``phi_cylinder_vacuum`` and ``phi_cylinder``, in
``tests/gates/derivations/_characteristic_point_mp.py``, are the routes the
characteristic reading and lift gates compare a cylinder's point flux
with: a solid cylinder under vacuum, and under a specular law, read at a
radius, written in mpmath from the closed forms. They are the first
clients outside ``orpheus/``, and the first whose value is a test's
oracle rather than a reference's reading. Each carries ``@traced_memo``
alone.

*Why these two, and no other route of the module.* `[M]` 2026-10-17, a
cProfile of a warm run of the slow rows
(``scratch/characteristic_architecture/p3/c3/drop/timing/profile/warm_slow.prof``):
82 % of the run was in ``mpmath.quad``, and 993 s of it in these two
routes, ``phi_cylinder_vacuum`` at 12 distinct points a run, 43.6 s each,
and ``phi_cylinder`` at 5, 94 s each. Each sat under ``lru_cache``, which
lives for one process, so every run paid for them again. The module's
other routes took milliseconds a point and kept their ``lru_cache``, apart
from the two stored after them (the next paragraph). That
split is the rule the exact medium's section set: a stored hit costs
about 8 ms in process (`[M]` the elegance review's probe3: 17.7 ms on the
first lookup, then 8.3 and 7.5 ms) and a first call a miss of at least
one interpreter start, about 0.5 s (:ref:`verification-reference-cache-miss-cost`),
so a millisecond route would cost about its own time on every lookup and
far more on its miss.

*The two routes stored after them* (``76d0409d``, the user's ruling of
2026-10-17, "store them"). ``specular_cylinder_total``, the block total of
a solid cylinder behind a partial mirror, now carries ``@traced_memo`` in
place of its ``lru_cache``, two keys a run. ``cylinder_polar_density`` is
new: the lift's polar-column rows computed it inline, the density
:math:`g(\theta) = \tfrac{1}{\pi}\sin\theta\int\psi\,d\alpha` at each
polar node of the reading's rule. It is stored as one column per call, not
one entry per node: a node costs 0.07 to 0.21 s, less than a miss, while a
column costs 9 to 27 s. `[M]` against the inline computation under
``bypass()``, the 12 nodes of the interface row: 9 bit-equal and 3 off by
1 ulp, against the row's band of 1e-13 (its tangency breaks are now
computed at 20 digits rather than 15). Over the 27 rows that read them
(``scratch/characteristic_architecture/p3/c3/concept2/rows_cold.log`` and
``rows_warm.log``, one run each, serially, under ``-O``), a cold store
took 162 s and a warm one 29 s; the polar row
``test_each_polar_node_of_the_reading_is_its_mpmath_density[wall-mirror]``
went from 14.3 s warm on the tree before (the scoping's measurement) to
1.5 s, and the partial-mirror row ``[0.99]`` of
``test_on_a_cylinders_partial_mirror_the_block_and_the_wall_reading_are_the_closed_forms``
read 83 s cold and 2.7 s warm. The census gate below reads the whole
module, so it covers both.

*Why this is not the size threshold refused for the transport block.* On
2026-10-15 the user ruled that ``group_transport`` is memoised whatever its
size (:ref:`characteristic-block-memo`): a threshold would have sent the
fast gates' cheap blocks through an in-process arm while production's
costly blocks took the memo arm, two routes to one value, chosen at the
call site by a performance number. Here no value has two routes. Each
oracle value has one route, its mpmath body, stored or not, and nothing in
production reaches it. Choosing which routes to store chooses which
functions are clients; it never chooses, call by call, which arm of one
function a value takes.

*Why no* ``lru_cache`` *sits above the memo.* A function cache stacked above
a memo serves, for the life of the process, whatever the memo returned
first, so under :func:`~orpheus.numerics.traced_memo.bypass` it serves the
value from before the bypass. `[M]` the elegance review's probe2: after a
module constant was patched under ``bypass()``, a memo under ``lru_cache``
returned 2.0 before and 2.0 after the patch, a bare memo 2.0 and 10.0.
The bypass is the spelling for a test that patches what a generation runs
(:ref:`verification-reference-cache-bypass`), so the stacking breaks it
silently and leaks a value between tests in either direction. At about
8 ms a hit the function cache buys nothing, so it was removed. The first
version kept it, above the memo, and exposed the child-resolution defect
of the process boundary (the section "The process boundary" above).

*The oracle's independence, gated.* The oracle is evidence about the
characteristic reference only if it shares nothing with it above the
trusted-library line (``instrument-doctrine`` X4), and the memo is the
oracle's first ``orpheus`` import. The gate
``test_each_stored_reference_imports_nothing_but_the_memo`` (foundation, in
``tests/gates/derivations/test_characteristic_reading.py``) is an AST
census, one row per stored independent reference: ``[mpmath-routes]``
and, since concept 3, ``[ps1982]``, the PS-1982 solver's module. ``tests.gates._corroboration.references_in`` lists every name a
module imports at any depth of its source, a function body's import
included; the gate reads the oracle module and every ``tests`` module it
reaches through its imports, at any depth, and allows one ``orpheus``
name, ``orpheus.numerics.traced_memo`` (with its ``traced_memo``). Its
control, in the same gate: the same census of the characteristic
reference's ``lines`` module must find ``orpheus`` imports beyond the
memo, or the gate fails with "the census sees no import". Its first reds,
each on a copy restored with a clean diff:

* an import of ``orpheus.numerics.quadrature`` at the top of the oracle
  module;
* an import of ``characteristic.lines`` inside one of its functions;
* an import of ``orpheus.geometry.boundary`` in ``_characteristic_mp``,
  which the oracle imports.

*The refuted first version.* The first gate compared the ``orpheus``
modules that fresh interpreters hold after importing the memo and after
importing the oracle. It is refuted for the question "does the oracle
import the code under test?", and both reviews found it blind twice. (1)
Importing the memo alone loads 76 ``orpheus`` modules (`[M]` 2026-10-17,
a fresh interpreter: 54 of ``orpheus.numerics``, whose package
``__init__`` imports it eagerly; 21 of ``orpheus.geometry``, the package
``orpheus.geometry.boundary`` whole among them; and ``orpheus``), every one
inside the reference's own import closure, so an added import of any of
them changed nothing the comparison saw. (2) An import inside a function
does not run at import time, so no interpreter holds it. The fact the
refutation establishes: what an import adds to ``sys.modules`` is not
whether the source imports it, and only the source can answer the second.

*What an entry pins, and what an edit there costs.* A stored entry pins
what its generation ran, and the generation imports the memo, so the
entry pins the skeletons of those 76 modules (`[M]` qa's review: 76
``orpheus`` module pins and 39 def pins of ``orpheus`` functions run at
import, beside the distribution pins of numpy and scipy). An edit to an
import, a constant or a signature anywhere in ``orpheus.numerics``'s
loaded modules or in the geometry they load (``orpheus.geometry.boundary``
included) regenerates both routes, about 993 s at the profile's
per-point figures (`[R]`). It costs minutes and never changes a value: by
the census, the routes import nothing of those modules beyond the memo.
The root is #605 (importing the memo imports the whole numerics package);
its fix shrinks the cone to the memo's own modules.

*What it bought.* `[M]` 2026-10-17
(``scratch/characteristic_architecture/p3/c3/oracle/summary.txt``, driver
``run.sh`` beside it): the 23 slow rows (``-m slow``) of
``test_characteristic_lift.py`` and ``test_characteristic_reading.py``,
each tree run once, serially, under ``-O`` and ``OMP_NUM_THREADS=1``, in a
fresh worktree with an empty store, cold and then warm.

.. list-table:: The slow rows of the lift and reading gates, before and
   after the oracle routes became clients (single runs)
   :header-rows: 1
   :widths: 40 20 20 20

   * - Tree
     - Cold
     - Warm
     - Entries in the store
   * - ``b2037667``, the routes under ``lru_cache``
     - 1 075 s
     - 353 s
     - 17
   * - ``333268fb``, the routes stored
     - 1 044 s
     - 127 s
     - 34

The 17 added entries are the 12 and 5 distinct oracle points. A warm run
saves 226 s. A cold run computes each route once either way, so its
31 s difference is between two single runs and is not a measured gain.
On ``333268fb`` the rest of the warm 127 s was computed in process, not in
the oracle (the ``--durations`` of ``after_warm.log`` beside the summary):
the characteristic reading's point read through ``point_flux``, the top
rows 28.9 s
(``test_the_reading_at_a_general_point_of_a_cylinder_is_the_mpmath_route[vacuum-wall]``)
and 25.6 s
(``test_a_cylinders_reading_is_the_frozen_mpmath_route[F4-mr3-vacuum-1.1]``),
and the lift's polar-column quadrature, 14.3 s
(``test_each_polar_node_of_the_reading_is_its_mpmath_density[wall-mirror]``).
The polar column is a stored route since ``76d0409d`` (above), and the
point readings are served by the reading's memo since concept 2 (below).

**The slow cylinder point readings, through the door** (#405 P3 concept
2 of ``.claude/plans/verification_tree_refinement.md``, 2026-10-17,
``94143ecd``; the user's ruling of 2026-10-17, "Door, with batteries in
bypass"). Eight ``slow`` rows of the characteristic gates read their point
flux as a question put to the reference,
``CharacteristicDerivation.evaluate(PointValue(x, g))``, one call per
group, instead of ``GalerkinSystem.point_flux`` in the test process. They
are the three rows of
``test_the_reading_at_a_general_point_of_a_cylinder_is_the_mpmath_route``
(``[vacuum-1.1]``, ``[vacuum-wall]``, ``[partial-1.1]``),
``test_a_cylinders_interface_point_is_read_to_the_bar``, the wall leg of
the two rows of
``test_on_a_cylinders_partial_mirror_the_block_and_the_wall_reading_are_the_closed_forms``
(``[0.9]``, ``[0.99]``; the block-total leg is unchanged),
``test_a_small_cylinders_layer_is_graded_outside_slow[wall]`` in
``test_characteristic_reading.py``, and
``test_a_cylinders_reading_is_the_frozen_mpmath_route[F4-mr3-vacuum-1.1]``
in ``test_characteristic_lift.py``. Each row's claim, route side and band
are unchanged.

*The posing.* ``_given_by_the_door``, beside the in-process ``_given`` in
``tests/gates/derivations/test_characteristic_reading.py``, takes
``_given``'s arguments and returns its type. It poses a
``GeometrySpecification``: a solid ``StructuredGeometry`` behind one outer
tag (bound in the ``exterior`` to the vacuum wall at albedo 0, a partial specular wall otherwise), pure
absorbers of the row's total cross sections, and the question
``FixedSource(RegionwiseConstant(q))``. It reads that specification at the
``Resolution`` that ``_door_resolution`` maps the row's (degree, layers,
line points) to: the grading ratio 0.4 and the least source points
:math:`2(p + 1)`. With no scattering and no fission the emission is the
source, so the door reads the quantity ``_given`` reads. The helper
refuses with a ``ValueError`` a hollow body, more than one law and a law
with a diffuse part; the door itself poses a hollow body and a diffuse
wall, so the edge is the helper's, not the door's. A test picks its reader per parameter with ``read=``
(on ``_check_route``, and on the lift's ``_reading``), so a test that
has both kinds of row carries both: in the lift file one set,
``_F_BY_THE_DOOR``, decides both a row's tier and its reader. The door
rows declare ``rests_on`` the door's own gate on a pure absorber,
``test_a_pure_absorbers_point_values_are_the_mpmath_route`` (E4),
per parameter, so the in-process rows of the same tests do not.

*Why the door, and not a memo of* ``point_flux``. A reading of the
reference is already a client: ``evaluate`` is keyed on the derivation
(the specification and the resolution) and the observable, and its solve
is a child entry. A point reading posed as a question is therefore served
by machinery that exists, and concept 2 memoises nothing new.
``GalerkinSystem.point_flux`` cannot be a client: ``GalerkinSystem`` is a
``dataclass(frozen=True, eq=False)``, compared by identity, so it has no
content, and the content encoder refuses it with ``ContentlessError``
(``orpheus.values.content``; the scoping's census). A test-side memo of
``_given`` would encode (its arguments are tuples, a string and a float),
but it would be a second spelling of the door, and it would put a hit or a
miss on the fast sphere and slab rows, which take milliseconds in process.
So the rows were re-posed as questions rather than given a memo of their
own. The open follow-up is one posing for both readers (#609): the system
``_given`` builds and the specification ``_given_by_the_door`` poses are
two spellings of one problem, which agree on the admitted domain.

*One ulp.* `[M]` the test-architect (the concept-2 specification's
section 2, ``ta/equiv_*.log``, on ``0da80c8c``, under ``bypass()``, each
row in a fresh process): 13 of 13 group readings of the eight rows are
within 1 ulp of ``_given``'s, 7 at 1 ulp and 6 equal, the largest
relative difference :math:`2.2 \times 10^{-16}`. ``orpheus/`` is the same
on ``0da80c8c`` and ``94143ecd`` (``git diff --stat``), so the
measurement carries to the landed tree. The ulp is the rounding of the
source's least solution (``source_emission``), which on a pure absorber is
the source in exact arithmetic: qa measured on
``test_a_small_cylinders_layer_is_graded_outside_slow[wall]`` (1 row of
the 8) that the source coefficients and the restriction are bit-equal and
that ``source_emission`` moves 5 of 12 emission entries by 1 ulp. By the
triangle inequality a row's error to its route moves by at most
:math:`2.2 \times 10^{-16}`, against measured errors of
:math:`2.4 \times 10^{-15}` to :math:`2.8 \times 10^{-14}` below a band of
:math:`10^{-12}`; the frozen lift row reads :math:`-4.22 \times 10^{-15}`
through the door and :math:`-4.33 \times 10^{-15}` through ``_given``,
against :math:`10^{-13}`. The equivalence is not the evidence: each row is
compared with its mpmath route, which shares nothing with the reference
above the trusted-library line.

*What a warm store no longer exercises in the test process.* On a warm
store the reading code these rows reach (``PointRule.of``,
``PointRule.row``, ``AttenuationRule``, ``GalerkinSystem.group_point_flux``)
does not run in the test process: the hit serves what that source
computed, as the block memo already does for the blocks the rows read.
Three consequences follow.

#. **Coverage.** A coverage capture taken on a warm store cannot
   corroborate a ``verifies`` claim from these rows: it reads ``refuted``
   or ``out_of_capture``. ``characteristic-quadrature`` keeps verifiers
   read in process. ``peierls-greens-cylinder-architecture`` does not: its
   only verifier is
   ``test_the_reading_at_a_general_point_of_a_cylinder_is_the_mpmath_route``,
   a door row, so its claim is corroborated only by a capture on a cold
   store or under ``bypass()`` (qa's review of concept 2, finding 1).
#. **Mutation.** A battery or a ``catches`` re-drop reaches these rows
   only under ``bypass()`` or on a mutated copy of the tree (the bypass,
   above). The concept-2 battery
   (``scratch/characteristic_architecture/p3/c3/concept2/ta/concept2_arms.py``)
   enters ``bypass()`` at ``pytest_configure`` and lists as ``UNREACHED``
   every row that ran under any other target. `[M]` the test-architect,
   on ``0da80c8c`` with the draft, over 13 rows (the 8 door rows, the
   small cylinder's interface row, the lift's F1 and F1b, and E4's two
   sphere rows), ``UNREACHED`` empty on every arm:

   .. list-table:: The concept-2 battery under ``bypass()``
      :header-rows: 1
      :widths: 34 33 33

      * - Arm
        - Red
        - Green
      * - none
        - none
        - 13 of 13
      * - the polar sine dropped from the point rule (the positive
          control)
        - 11 of 13, at 0.36 to 0.90
        - E4's two sphere rows (a sphere has no polar sine)
      * - the door's answer scaled by :math:`1 + 10^{-9}` (the route
          witness)
        - 10 of 13, at :math:`10^{-9}`: the 8 door rows and E4's two
        - the small cylinder's interface row, F1 and F1b, exactly the
          rows read in process
      * - the layer term removed (ERR-105)
        - the small cylinder's wall :math:`4.8 \times 10^{-7}` and
          interface :math:`1.1 \times 10^{-8}`
        - —

   The in-process default tier lost no first red. `[M]` qa, on the landed
   tree without ``bypass()``, ``-m "not slow"``: the arm grading the
   layer only for :math:`0 < a < 1`, whose only red in the battery was a
   door row, reddens the small cylinder's interface row and the lift's
   ``F2-thick-wall`` and ``F3-thin-wall``; the arm grading the outermost
   panel only, likewise, reddens the small cylinder's interface row; the
   arm none, 6 of 6 green. One documented first red is blind on today's
   tree whichever reader the row uses: the level not passed by the line
   set, on the partial-mirror wall rows, leaves both green under
   ``bypass()`` (#608; ``_given`` reads the same quantity within 1 ulp).
#. **What the door rows exercise instead.** On a cold store or under
   ``bypass()`` they also run the door's posing: the walls of a law
   (``Walls.of``), the coefficients of a ``RegionwiseConstant``,
   ``source_emission``, ``solve`` and the ``PointValue`` arm of
   ``evaluate``.

*Why one row stays in process.*
``test_the_volume_integral_of_the_reading_is_the_galerkin_block`` (C6)
integrates the reading OPERATOR, one row per point of a 648-point volume
rule against every basis column, against the block. The door poses a
source as a ``RegionwiseConstant``, which cannot pose one basis function,
so C6 has no question to put through it, and it reads in process. Its two
``slow`` white rows are the slowest of the warm ``slow`` run below.

*What it bought.* `[M]` 2026-10-17
(``scratch/characteristic_architecture/p3/c3/concept2/landed_cold.log``,
``landed_warm.log`` and ``landed_warm_slow.log``): the landed tree,
``test_characteristic_lift.py`` and ``test_characteristic_reading.py``,
on the main tree's store, each run once, serially, under ``-O`` and
``OMP_NUM_THREADS=1``.

.. list-table:: The lift and reading gates after concept 2 (single runs,
   the main tree's store)
   :header-rows: 1
   :widths: 40 30 30

   * - Selection
     - Store
     - Time
   * - both files, 240 rows
     - cold for the door rows' keys
     - 1 882 s
   * - both files, 240 rows
     - warm
     - 217 s
   * - the 23 ``slow`` rows
     - warm
     - 25.8 s

Warm, each door row takes 0.03 to 0.13 s; the slowest ``slow`` rows are
C6's two white rows, 10.9 s each, and the polar-column row
``[wall-mirror]``, 1.5 s. This is not the protocol of the concept-1 table
above (a fresh worktree with an empty store), so no before and after ratio
against its 127 s is claimed. 240 of 240 collected ids are unchanged, and
217 of 217 under ``-m "not slow"`` (qa, ``--collect-only``).

**The PS-1982 vacuum-sphere solver** (#405 P3 concept 3, 2026-10-17;
the user's ruling of 2026-10-17, "Memo the solver").
``solve_ps1982_vacuum_sphere``, in
``orpheus/derivations/continuous/peierls_nystrom/ps1982_reference.py``, is
the independent reference the characteristic reference's one-group
vacuum sphere is compared with
(``tests/gates/derivations/test_characteristic_independent_references.py``:
:math:`k` to :math:`10^{-5}` relative on four cases, the flux shape to
:math:`10^{-4}` on the first two). It solves Pomraning and Siewert's
scalar-flux integral equation :cite:`PomraningSiewert1982` (Eq. (21)) by
power iteration on 40 Gauss–Legendre nodes, each application of the
operator an adaptive ``scipy.integrate.quad`` per node over the kernel
:math:`E_1(|r-x|) - E_1(r+x)` with a break at :math:`x = r`
(:doc:`/theory/references/peierls`, "The Pomraning-Siewert 1982 vacuum
sphere"). It carries ``@traced_memo`` alone.

*Why it.* `[M]` 2026-10-17
(``scratch/characteristic_architecture/p3/c3/concept3/scoping.md``, its
section 4): after concept 2 a warm ``-m slow`` run of the 40 slow rows in
their 10 files took 82.9 s, and 55.8 s of it was the PS-1982 rows
(21.74, 21.22, 12.11 and 0.73 s), all in the solver: a full run makes
four distinct solves (three under ``-m slow``; the shape rows read the
same cases as the first two :math:`k` rows), each 11.6 to 20.7 s at the
gates' 40 nodes, 59 or 110 iterations (`[M]` ``probe_calls.log``, one run
each, under ``-O`` and ``OMP_NUM_THREADS=1``). The solver is a client the
memo admits as it stands. It has one caller, the test helper ``_ps1982``,
and no production caller (the Nexus graph's transitive callers and
``git grep``). It is pure: no module state, a deterministic start, a
repeat call bit-equal in :math:`k` and the flux, and no warning recorded
on any of the four cases under ``simplefilter("always")`` (its control,
``quad(1/x, limit=5)``, records one), so #607 is inert here. Its result,
``PS1982Result``, is a frozen dataclass of seven plain-data fields (three
floats, two float arrays, an integer and a flag) that the return
annotation names, so the payload rebuilds it through its constructor. A
solve that exhausts its iterations is a value and would be stored; the
``converged`` assertion stays in the helper. The key is exact, so
``R=5`` and ``R=5.0`` are two entries; the helper passes floats.

*Why on the solver, and not on the test helper.* A memo on ``_ps1982``
would pin the test module's import closure, which is the whole
characteristic reference, so every edit to the code under test would
regenerate its independent reference; and it would be a second spelling
of the solver. That arm is refuted. On the solver the memo is one route
per value, the precedent of the multi-region solvers above. The helper's
``lru_cache`` was removed: a function cache above a memo serves, for the
life of the process, the value from before a ``bypass()`` (the oracle
routes' paragraph above says why).

*The alternative not taken: memoise the eigenpair.* The eigenproblem
depends on the optical radius :math:`R\Sigma_t` alone: the body's power
iteration reads neither :math:`\Sigma_s` nor :math:`\nu\Sigma_f`, and only
its last three lines (the critical albedo, the denominator and
:math:`k`) do. `[M]` the scoping (section 2): two cases with the same
optical radius and different :math:`\Sigma_s` return bit-equal
``lambda_max``, ``phi``, ``r_nodes_optical`` and ``iterations``, and two
optical radii differ (the control). So the four cases are two
eigenproblems and four arithmetic conversions. A memoised eigenpair keyed
on the optical radius, the node count and the iteration budget, with
:math:`k` converted in process, would halve the cold cost, to about 36 s
against about 72 s, and let the default tier's entry serve a slow row's
eigenproblem. It was not taken: it restructures the reference code (a
second function and a second result type) and keys on a computed float
product, while warm the two arms cost the same. The user chose the
solver unchanged.

*What an entry pins, and what an edit there costs.* A generation imports
``ps1982_reference``, so an entry pins the module's functions it runs
(``solve_ps1982_vacuum_sphere``, ``_apply_L_operator`` and its integrand)
and the skeletons of about 90 ``orpheus`` modules: the 76 the memo itself
loads (#605) and 14 more that the package chain from ``orpheus.derivations``
down to ``ps1982_reference`` loads, the module itself among them (`[M]`
qa's review, ``scratch/characteristic_architecture/p3/c3/concept3/qa/``:
90 loaded, the memo's 76 all inside them, none under ``characteristic``).
An edit to the module, or to an import, a constant or a signature in that
cone, regenerates all four entries, about 72 s (`[R]`, from the per-solve
times); the values never change under those edits. A docstring edit does
not regenerate them: a def is pinned by its normalised AST and a module by
its skeleton, both with docstrings stripped (gate M2.1). `[M]` 2026-10-17
(``scratch/characteristic_architecture/p3/c3/concept3/probe_docedit.log``
and ``probe_moddoc.log``, each edit made on a copy-protected file and
restored): after a docstring-only edit to ``solve_ps1982_vacuum_sphere``
the file ran warm, 9 passed in 3.15 s, and after a module-docstring edit
in 2.50 s; no entry's ``DataPin`` names a ``.py`` file
(``probe_manifest.py``). [REFUTED 2026-10-17] the review's first reading,
that a docstring edit regenerated all four entries (a 227 s run): it was
inferred from file times, and that run's cause was not recovered.

*What a warm store no longer exercises.* The solver does not run on a
hit, so a mutation of it reaches the PS-1982 rows only under
``bypass()`` or on a mutated copy of the tree, as for the door rows
above. `[M]` 2026-10-17
(``scratch/characteristic_architecture/p3/c3/concept3/probe_bypass.py``,
re-run by the archivist): with ``_apply_L_operator`` scaled by 1.001, the
case ``tauR2.5-c0.4`` read :math:`k = 0.0563821101969148` from the store,
the unmutated value, and :math:`0.056463960765175535` under ``bypass()``,
a relative move of :math:`1.45 \times 10^{-3}`, 145 times the
:math:`10^{-5}` band. The PS-1982 rows carry no ``verifies`` or
``catches`` marker, so a coverage capture on a warm store refutes no
claim of theirs.

*What it bought.* `[M]` 2026-10-17
(``scratch/characteristic_architecture/p3/c3/concept3/landed_cold.log``
and ``landed_warm3_1.log``, the warm run re-captured after the review's
edits): the nine rows of
``test_characteristic_independent_references.py``, the slow ones
included, on the main tree's store, each run once, serially, under
``-O`` and ``OMP_NUM_THREADS=1``: 87 s with the store cold for the four
PS-1982 keys, 2.29 s warm (2.25 s on a repeat), each PS-1982 :math:`k`
row 0.02 to 0.22 s warm. The 40 slow rows of the 10 characteristic files
(``slow10_first.log`` and ``slow10_warm.log`` beside them, the same
protocol): 216 s on the pass that regenerated the entries the review's
edits had invalidated, then 27.1 s warm, against 82.9 s before this
concept; 20.7 s of the 27.1 s is C6's two white rows (10.35 and 10.32 s),
which read in process.

*The 24 files re-timed after concepts 1 to 3.* `[M]` 2026-10-17
(``scratch/characteristic_architecture/p3/close/timing/``: ``run.sh``,
``results.tsv``, ``logs/`` with ``--durations=0``). The protocol is the
one of 2026-10-14: the same 24 files, one at a time, under ``-O`` and
``OMP_NUM_THREADS=1``, in a detached worktree of ``main`` at
``a456a9bf``, in four passes. The passes are: the memo bypassed; the
store emptied (cold); the store warm (warm1); and warm with
``-m "not slow"`` (warm2). 96 of 96 rows pass. The 3 warm2 rows that exit
5 are files whose every test is ``slow``.

.. list-table::
   :header-rows: 1
   :widths: 40 15 15 15 15

   * - file (seconds)
     - bypass
     - cold
     - warm1, 2026-10-14
     - warm1, now
   * - ``test_characteristic_assembly``
     - 2 609
     - 2 622
     - 3 535
     - 2 579
   * - ``test_characteristic_system``
     - 845
     - 892
     - 1 105
     - 4
   * - ``test_characteristic_lift``
     - 642
     - 638
     - 1 090
     - 134
   * - ``test_characteristic_albedos``
     - 744
     - 1 078
     - 988
     - 6
   * - ``test_characteristic_reading``
     - 870
     - 1 170
     - 1 032
     - 88
   * - ``test_characteristic_convergence``
     - 712
     - 819
     - 933
     - 3
   * - ``test_characteristic_transport``
     - 136
     - 136
     - 166
     - 138
   * - ``test_characteristic_independent_references``
     - 210
     - 223
     - 213
     - 4
   * - ``test_characteristic_reference``
     - 83
     - 244
     - 107
     - 15
   * - all 24 files
     - 8 630
     - 8 571
     - 9 307
     - 3 073

warm2 took 487 s, against 991 s on 2026-10-14. The bypass and cold
passes fell from 11 601 s and 10 329 s; `[R]` through the reference's own
speed-ups since (the axial lift, the derived node counts), not measured
per cause.

What remains warm is the system under test, computed in process by
design. ``test_characteristic_assembly`` holds 2 579 of the 3 073 s
warm. Of that, 2 493 s is its 21 ``slow`` rows (warm1 minus warm2), and
12 cylinder rows account for 2 345 s. The 6 closed-body conservation rows
on the three-region cylinders take 277 to 418 s each, and the 5 cylinder
escape and transmission closed forms take 30 to 97 s each, and the
interface between equal materials on a cylinder takes 109 s. These gates
build the transport block through ``LineRule.of(...).transport``. They do
not go through the stored ``group_transport``, because the block is what
they verify, and a served block would make them read the store. So no
stored reference can make them cheaper. Two things can: a cheaper
cylinder block (#587, the non-tensor (b, θ) line rule), or re-running a
verdict only when the code it rests on changes, which is ruling 15 of
the verification-tree plan, the same route as C6. The remaining warm
cost of ``lift``, ``transport`` and ``reading`` (88 to 138 s each, of
which 63 to 137 s is outside ``slow``) is likewise in-process transport:
the Volterra-triangle rows lead the durations of ``lift`` and
``transport``, and the fundamentals' readings and C6 those of ``reading``.

**A solve that does not converge.** The solver's result is a value
whatever its convergence, so an exhausted solve is a cached value; the
reading of it is a refusal, and a refusal is an exception, which writes no
entry (M4.8). The same question is open for warnings (#607): the caller's
warning filters do not reach a generation, which runs in another process
under the interpreter's default filters, so a warning that a test module's
``filterwarnings("error::RuntimeWarning")`` would raise in process is only
printed in the generating child, and the entry is stored. `[R]` inert
today: no stored generator's correctness rests on a warning raised as an
error (`[M]` 2026-10-17, qa's grep of their module tree). Whether the
filters become part of what a generation inherits, and so of the key, or
the page states that an entry is a value whatever the caller's filters,
is #607's ruling.

.. _verification-reference-cache-miss-cost:

What a miss costs
-----------------

`[M]` 2026-10-15 (#602; the explorer's profile and its probes under
``scratch/characteristic_architecture/p3/m602/``). The fixture is five
cheap ``group_transport`` blocks (three one-region and two three-region
spheres, degree 3), each repeat salted so that every call misses, each
repeat on a temporary store, ``-O``, serial. Before #602 a miss was 2.65
to 2.88 s, of which:

.. list-table:: One cheap miss before #602 (median, ms)
   :header-rows: 1
   :widths: 60 20

   * - Part
     - Time
   * - interpreter start, to the boot script's first line
     - 37
   * - imports (``orpheus.numerics``, then the client's package)
     - 467
   * - the computation itself (one-region sphere)
     - 30
   * - the manifest
     - 1 791
   * - the child's teardown after its answer is written
     - 369
   * - the parent's read-back check (about 800 on a process's first)
     - 14

The recorder itself costs almost nothing: it is a ``sys.monitoring``
``PY_START`` hook that disables itself after each code object's first
report, and it added 37 ms to the imports and nothing measurable to the
computation. The manifest and the read-back held four costs, each a fact
about unchanged files recomputed by every fresh process; the two rounds
of #602 removed them.

* **Every source file the run touched was parsed, and every def in it
  digested, in every miss.** Each of the 129 source files the run touched
  was parsed twice and stripped of docstrings twice, and all 1 948 of
  their defs were dumped and digested, when 247 ran. Round one
  (``8e02236f``, ``e666e58d``) parsed each file once per process and
  digested a def only when a pin asked for it. Round two (``155fbdbd``)
  keeps the answer on disk: a text's index is built from one parse by the
  first generation that meets it and read back by every later process, so
  a miss parses only text no generation has indexed, and a process's first
  validation reads indexes instead of parsing every pinned file
  (:ref:`verification-reference-cache-source-index`). An index on disk
  cannot know which defs later pins will ask for, so it digests every def
  when it is built, once per content. The skeleton is the digest of a
  shallow copy of the module's statement spine in which every outermost
  def's body reads as ``pass``, so the parsed tree is never changed; the
  docstring strip descends statements only (a docstring is the first
  statement of a module, class or def body), where it walked every
  expression node.
* **A distribution's version was looked up per code object**: 1 245
  lookups for 6 distributions, 0.64 s, each re-reading the metadata.
  Round one looked it up once per distribution and per checkout (one
  manifest or one validation); since 2026-10-17 it is looked up by
  ``_Installed.installed`` (``_Installed.version`` until cut E, which added
  the ``RECORD`` digest) from ``_DISTRIBUTIONS``, a read kept for the
  process while what the find looked at is unchanged on disk
  (:ref:`verification-reference-cache-kept-reads`), so the manifest and
  the validation read a version by one definition. A distribution that owns a file the run touched and has no
  installed version is ``Unpinnable`` (``_distribution_pin``).
* **The owner of each site-packages file came from a top-level map**,
  ``importlib.metadata.packages_distributions``, which read every installed
  distribution's metadata, about 0.22 s per miss. The per-file table that
  replaced it for correctness (``0c569df8``) reads the ``RECORD`` files in
  18 to 20 ms.
* **The child freed its objects before exiting** (0.33 s). The boot
  script ends with ``os._exit(0)`` once the answer is written and both
  standard streams are flushed.

None of these changes what an entry pins. Over every tracked Python file
(1 234 files, 15 515 defs) round one's parse and the original give
identical skeleton and def digests, with controls that a body edit moves
only the def's digest and a constant edit only the skeleton's; the same
comparison after ``e666e58d`` (1 234 files, 15 522 defs) found 0
mismatches. After round two, ``_def_pin`` read from the index equals the
original tree-based ``_def_pin`` on 33 750 code objects of 1 236 files,
built cold and read back from disk, with 0 mismatches
(:ref:`verification-reference-cache-source-index`).

.. list-table:: One cheap miss, per round (seconds; each process's first miss excluded)
   :header-rows: 1
   :widths: 34 22 22 22

   * - Tree
     - One-region block
     - Three-region block
     - A process's first read-back validation
   * - Before #602
     - 2.65 to 2.88 (all five blocks)
     - (not measured apart)
     - about 0.8
   * - Round one, ``e666e58d``
     - 1.17 to 1.20
     - 1.33 to 1.35
     - 0.27
   * - Round two, ``155fbdbd``
     - 0.66 to 0.70
     - 0.84 to 0.85
     - 0.036

`[M]` 2026-10-15 for the last two rows, the logs
``run_after_loader.log`` and ``run_after_R2.log`` in the same directory,
five blocks times three repeats in one process, a warm bytecode cache and
(round two) warm source indexes, re-read by this chapter's probe
(``archivist_numbers.py`` and its log in the same directory, 2026-10-16). Each cell
is the range over every draw but the process's first; the per-block medians
are 1.18 to 1.20 and 1.34 to 1.35 s (round one) and 0.69 to 0.70 and
0.84 to 0.85 s (round two). A hit is 0.010 to 0.013 s in both rounds, and
a later read-back validation 13 to 39 ms (round one) and 11 to 13 ms
(round two). The blocks themselves take 0.02 to 0.17 s in process.

**The three files.** The 431 not-slow rows of the system, reading and
assembly gates, in a fresh worktree of each commit (an empty store), one
serial run each (`[M]` 2026-10-15,
``scratch/characteristic_architecture/p3/m602/timing/summary.txt`` and
its logs, ``r2_*.log`` for round two):

.. list-table:: The system, reading and assembly gates after #602
   :header-rows: 1
   :widths: 28 30 42

   * - Run
     - Round one, ``8e02236f``
     - Round two, ``155fbdbd``
   * - No memo (``bypass()`` entered for the session)
     - 234 s, 244 s (two runs)
     - 236 s
   * - Cold store
     - 538 s; 257 block entries, 20 MB
     - 406 s; 257 block entries (19.9 MB), 129 source indexes (0.4 MB),
       127 bytecode entries (3.2 MB)
   * - Warm store
     - 177 s, 177 s (two runs); no entry written
     - 168 s; no entry written

Before #602 the cold pass was 1 303 s, on the block memo's prototype,
which stored 383 entries (:ref:`characteristic-block-memo`). The cold pass
now costs 170 s more than no memo, about 0.66 s per entry, where round one
cost about 300 s, about 1.2 s per entry. The warm pass saves 68 s of
236 s, less than the prototype's 102 s warm pass suggested, and the
per-test durations say why. The landed assembly gates build their blocks
in process through ``LineRule.transport`` (the slowest, a cylinder block
under a partial mirror, takes 16.8 s warm and without the memo alike), and
the reading gates' point rows and their ``bypass()`` rows are also
computed. The prototype had routed the assembly gates' blocks through the
memo, which is why it stored 383 entries and ran warm in 102 s.

The three files were re-timed on 2026-10-17, at ``main`` ``3ff7866b``
(221 s without the memo, 816 s cold, 222 s warm) and on the first version
of the fix that keeps a validation's reads for the process (227 s, 725 s,
155 s), and on the landed fix, ``8fbb1d61`` (221 s, 820 s, 163 s); the
table and its reading are in
:ref:`verification-reference-cache-kept-reads`.

**What the two rounds left.** After round two a cheap miss is the floor of
a fresh process plus the computation and a manifest built from indexes:
the interpreter start and the imports, about 0.5 s (37 ms and 467 ms in
the first table), are paid again by every miss, whatever the manifest
costs. Only a process that has already imported the package can remove
them, and that is cut E, the fork template, below.


.. _verification-reference-cache-fork-template:

The fork template (#602, cut E)
-------------------------------

Where forking has been measured safe (the scope boundary below), a miss is
no longer generated by a new interpreter. It is generated by a child
forked from a **template**: an interpreter the asking process starts once,
which starts exactly as a generating process does, does the imports, and
then forks one generating child per miss. The child begins with the
imports already done, so the interpreter's start and the imports are
paid once per template instead of once per miss. The design is cut E of
the round-two specification
(``scratch/characteristic_architecture/p3/m602/round2_gate_spec.md`` §3);
the code is :func:`~orpheus.numerics._traced_memo_boot.template` and, in
the asking process, ``_SLOT.answer`` and ``_Template`` in
:mod:`orpheus.numerics.traced_memo`; the gates are
``tests/gates/numerics/test_traced_memo_template.py``.

**Why a forked child may stand in for a fresh process.** The ruling of
:ref:`verification-reference-cache-fresh-process` is about what the
generating process holds in memory, not about how it was started. A
child forked from the template holds what the template held at the fork,
so the template is built to hold nothing a fresh process would not:

* **The template is a fresh interpreter, never a fork of the asking
  process.** It is started with ``subprocess.Popen`` by the same command
  line as a fresh generating process, so no object, cache or monkeypatch
  of the asking process reaches it (hazard 3; gate E4 patches the
  generator's helper in the asking process before the template starts,
  and the forked entry holds the honest answer).
* **The template runs no generation.** It only imports; each job is run
  by a child, after the fork, and the child exits when it has answered. A
  ``functools.cache`` warmed by one generation therefore dies with that
  child and is never inherited by the next (hazards 1 and 2; gate E2: a
  cached helper run by two consecutive forked generations is pinned by
  both). The red on record is battery arm E-b, a template that runs each
  job itself before forking it, which reddens E1, E2 and 29 other rows
  (31 of 196 in battery F).
* **The child's recording is the fresh process's.** The template starts
  recording from its first line, exactly as a generating process does
  (the shared start, ``_boot``: the recording, the pinned-source loader
  of :ref:`verification-reference-cache-pinned-source-loader`, the
  reserved answer descriptor, the envelope carrying ``sys.path``, the
  memo's import). When the imports are done it takes a **snapshot** of
  the recording (``Recording.snapshot``). Each forked child puts the
  snapshot back (``Recording.restore``) and re-enables the
  ``sys.monitoring`` events of every code object
  (``sys.monitoring.restart_events``: a code object reports once and then
  disables itself, so without the re-enable a function the template had
  already started would never report in the child). The child's recording
  then holds what a fresh process's held after the same imports, and the
  generation adds to it from there. The fields copied are one tuple,
  ``_RUN_STATE``, so a field is carried by the snapshot and the restore or
  by neither. Gate E1 generates three misses forked from one template and
  the same three in fresh interpreters: the pin sets agree entry by entry
  (less the memo's own def pins, which the memo pin covers and which
  differ because the template's job loop runs memo code), the payload
  digests agree, and the three forked entries pin the same code.
  `[M]` 2026-10-16, battery E arm E-a: a child that does not restore the
  snapshot pins the same set on E1's generator, because everything the
  template runs after the snapshot is memo code or the standard library;
  the restore keeps the memo's own pins exact, and no gate reddens on its
  removal.
* **The child cannot tell the two routes apart.** ``sys.argv`` is cut to
  the script alone before anything is imported (the template is started
  with ``--template``; a generator reading ``len(sys.argv)`` saw 2 forked
  and 1 fresh, qa's review of cut E, finding 3), and the child's standard
  input is ``/dev/null``, at its end at once, as a fresh process's input is
  once its job has been read (a child that kept the template's job pipe
  as its input blocked on it, finding 4). Gate ``test_qa3_qa4_…`` asserts
  both on both routes; battery F arms R-3 and R-4, each removing one of
  the two lines, redden it alone.

**The template re-validates its imports before every fork.** A template
lives across many misses, and the files its imports read can change while
it lives. Its module state was built from the files as they were at its
start; a child forked after a change would run that old state while its
manifest, written after the generation, pins the files as they are now:
a stale answer under current pins, an entry that validates and is wrong.
So the template pins its own imports as a generation pins its run: when
the imports are done it builds the manifest of the snapshot
(``_manifest``), adds a data pin on the whole bytes of every first-party
source the loader compiled (``Manifest.with_compiled_sources``), and
before each fork it validates that manifest
(:func:`~orpheus.numerics.traced_memo.validate`). The compiled sources are
pinned whole, not only by the def and skeleton pins of the imports' own
run, because that run pins a def only if its body ran: the body of a def
the imports defined and never called has no pin in it, yet a generation
forked from the template later calls that def and runs its body as the
template compiled it. What the check
covers, each with its gate:

* **An edit to a module the template imported** (gate E3, two rows: the
  body of a def the import never ran, and a module constant). The next
  miss returns the new code's answer from a new template, and its entry
  validates. The red: a template that does not re-check the bytes it
  compiled forks the old code under the new pins.
* **A data file, a directory listing or a distribution the imports
  read** (gate E9: a module the template imports reads a first-party
  data file at import into a module constant; after one miss the file is
  edited; the next miss returns the new value from a new template, and
  the refusal names the data file). `[M]` 2026-10-16, battery F arm E-n,
  the readiness check without the validation: 3 rows red, E9 and the two
  E3 rows.

**And it checks that it runs one OS thread.** A fork copies the calling
thread only, together with every lock the other threads held at that
moment; a lock held by a thread that does not exist in the child is held
for ever there. A BLAS library starts a pool of native threads at its
first call, and a fork after it crashes the child: `[M]` 2026-10-16, the
fork-safety protocol's positive control, a template that multiplies two
800 × 800 matrices before forking, 20 of 20 children killed by a signal
(on the specification's earlier probe, ``ta/p7_fork_safety.py``, 5 of 5
by ``SIGSEGV`` with ``VECLIB_MAXIMUM_THREADS`` unset). ``threading``
cannot see native threads (``threading.active_count`` reads 1 after a
BLAS call), so the count is the operating system's:
``_libproc_thread_counter`` reads ``threadnum`` from libproc's
``proc_pidinfo(PROC_PIDTASKINFO)``, and ``_fork_unsafe`` refuses any
count above 1. Gate E5 reads the counter at 1 in a fresh process, 2 with
one started ``threading.Thread``, and above 1 after a 1500 × 1500 matrix
product; battery F arm E-e, the refusal removed, reddens E5 and S4 (below).

**The answers are a closed sum.** Every frame a generating process or a
template sends back unpickles to one of four types in
:mod:`orpheus.numerics.traced_memo`:

* ``_Generated``: the entry is written;
* ``_Failed(error, described)``: the call, or the writing of its entry,
  raised; the exception crosses pickled, and described for the case the
  parent cannot rebuild it (as before E);
* ``_Refused(reasons)``: a template will not fork (its imports' manifest
  no longer validates, or it holds more than one OS thread); it exits
  after answering;
* ``_Died(status)``: the child a template forked ended without answering,
  with that wait status. The asking process raises ``RuntimeError``
  naming the death; the template lives on (gate E7c: a generator that
  calls ``os.abort()``; the next miss is answered by the same template).

The first two come from either route, the last two only from a template.
Every pipe carries length-prefixed frames (``frame``, ``_read_frame``,
``_write_frame``), so the template's one job pipe carries a sequence of
jobs. One value, ``_Launch``, holds how a generating process is started
on either route: the interpreter, the ``-O`` level, ``sys.path``, the
working directory and the environment handed to the job. A fresh process
is started from it; so is a template; and a template serves a job only
when its launch equals the job's, because the snapshot it forks was taken
under that state. A changed environment variable, working directory or
``sys.path`` therefore closes the template and starts another (gate E6:
each of the three changes between misses, and each next miss sees the
new state and validates).

**The template's lifetime.**

* **Replaced after 256 forks** (``_FORKS_PER_TEMPLATE``): the template's
  resident size grew about 16 kB per fork in the protocol (99 to 156 MB
  over 3 600 forks, `[M]` 2026-10-16). Gate E7d: 260 misses leave no
  zombie child and no growth of the template's open descriptors over its
  first 200 forks, and the template is replaced.
* **A refused or vanished template is replaced at the same miss.** A
  template that refuses after it has forked, or whose pipe has closed,
  is closed and a new one is started for the same job; the new one is
  fresh, so a refusal from it turns forking off (next item), and a miss
  that two templates leave unanswered is generated by a fresh process.
  Gate E7b kills the template between two misses; gates E3 and E9 see the
  first template refuse after an edit and the second answer.
* **A refusal by a FRESH template turns forking off for the process.**
  When the template's first job is refused, the cause is in the imports
  themselves (an OS thread an import starts, say), and it would hold for
  every template the process started. The reasons are kept as
  ``_SLOT.refused``, and every later miss of the process spawns a fresh
  interpreter. A child forked from the asking process keeps the refusal,
  because it inherits the imports that caused it. Gate S4: a module the template imports starts a thread;
  three misses are answered by ONE template start, not three; battery F
  arm R-S4, the refusal not remembered, reddens it alone. The choice is
  conservative, and the case it closes is concrete: an import that probes
  a file, finds it absent and creates it would make every new template's
  first validation fail (the presence pin reads ``absent``, the file is
  there), so templates would be started and refused at every miss, each
  start costing more than the fresh process it replaced. Turning forking
  off costs speed only: the fresh route is the route that ran before E,
  and an answer is never served from a template that refused.

  The refusal is warned once, when forking turns off, as
  :class:`~orpheus.numerics.traced_memo.TemplateRefused` with the
  template's reasons. Until then it was silent. A process that turned
  forking off ran slower, and its forked-route rows skipped, with no
  cause printed. ``[M]`` 2026-10-16: one full ``-m "not slow"`` run of the
  bypass drop's branch showed three more skips than ``main``; a re-run of
  the same commit matched ``main``; the cause was not recovered. Gate S4
  asserts the warning, and a mutation removing the warning reddens it.
* **A template ends with its parent.** It reads jobs until its job pipe
  ends, and a parent that exits, even by ``SIGKILL`` with no ``atexit``,
  closes that pipe (gate E7a: the template is gone within 5 s); a parent
  that exits normally closes it in an ``atexit`` handler.
* **The template state is one object.** ``_SLOT``, a ``_TemplateSlot``,
  holds the template, the lock below and the refusal above, with the
  verbs ``answer`` (start, replace or refuse as the launch, the fork count
  or the refusal requires, then return the child's answer), ``close``
  (at exit) and ``forget`` (in a forked child).
* **One miss at a time uses the pipes.** A lock (``_SLOT.lock``)
  serialises the asking process's use of its template: without it, four
  threads missing at once interleaved their frames and read each other's
  answers (`[M]` 2026-10-16, qa's review of cut E: 14 of 40 calls raised).
  Gate ``test_qa1_…``, twenty rounds of four threads released at once by
  a barrier. Battery F's arm R-1, the lock removed, left the row green
  (`[M]` 2026-10-16, ``ta/batteryF/verdicts.txt``); the row was revised
  after that run so that its generator sleeps 0.05 s and the four misses
  overlap in the pipes, and a red of the revised row under R-1 is not in
  the battery tables.
  An ``os.register_at_fork`` handler, ``_SLOT.forget``, runs in a child
  of the asking process: it closes the child's copies of the template's
  pipes (the template sees the end of its input only when every holder
  of its standard input has closed it, so a child that kept a copy would
  keep the template alive after the parent closed its own) and starts
  with no template and a new lock (a lock held by a thread the fork did
  not copy is never released in the child).
* **A nested miss spawns.** A generating process, whether fresh or
  forked, never forks its own nested misses: ``_forks_here`` is false
  wherever the pinned-source loader is installed, so a nested miss starts
  a fresh interpreter as it always has. Gate E8: an entry generated as a
  nested miss inside a forked generation pins exactly what the same
  call's direct entry pins.

**The scope boundary: the platform the protocol ran on.**
``FORK_MEASURED`` in :mod:`orpheus.numerics._traced_memo_boot` is one
table, the only door to forking: a platform, named by the operating
system, the machine and numpy's BLAS
(:func:`~orpheus.numerics.traced_memo.fork_platform`), forks only when the
table lists it, with that platform's OS-thread counter. Today it lists
one, ``("darwin", "arm64", "accelerate")``, the platform the fork-safety
protocol ran on. Every other platform, Linux with OpenBLAS on continuous
integration included, spawns a fresh interpreter per miss. It is tagged
``SCOPE-BOUNDARY[guard]`` in the code: the machinery that would widen it
is the protocol run on the other platform and an OS-thread counter for
it. The gate (``test_a_platform_the_fork_protocol_never_ran_on_spawns_fresh``)
reports an OpenBLAS platform and requires one fresh generation and no
template.

**The fork-safety protocol** (§3.3 of the specification; `[M]` 2026-10-16,
``scratch/characteristic_architecture/p3/m602/ta/fork/full/driver.log``,
driver ``ta/fork/driver.py``, 02:07 to 02:56, serial, ``-O``). A
measurement template imports what a generating process imports for the
characteristic clients (the memo, numpy, scipy and the characteristic
package), asserts one OS thread before each fork, and forks three
workloads: W1, a 1200 × 1200 matrix product, ``eigh``,
``scipy.linalg.solve``, ``spsolve``, ``quad`` and ``expn``; W2, the five
cheap ``group_transport`` blocks and one characteristic reading; W3, a
child that starts a thread pool itself. Each child's result digest is
compared with a fresh interpreter's under the same environment. The
matrix is ``OMP_NUM_THREADS`` × ``VECLIB_MAXIMUM_THREADS``, each unset, 1
or 8. The protocol measures the fork, not the production template's
other machinery: the recording, the snapshot, the validation and the
lifetime are the gates' (E1 to E9, above).

.. list-table:: The fork-safety protocol: nine cells, 3 600 forks each
   :header-rows: 1
   :widths: 40 60

   * - Quantity
     - Over the nine cells
   * - Forks answered and equal to a fresh interpreter's digest
     - W1 3 000 of 3 000, W2 300 of 300, W3 300 of 300 in every cell:
       32 400 of 32 400
   * - Signals, non-zero exits, hangs, digest mismatches, a second
       thread before a fork
     - 0 in every cell
   * - New Python crash reports (``~/Library/Logs/DiagnosticReports``)
     - none
   * - OS threads of the template at its start and after its imports
     - 1 and 1 in every cell
   * - A child's time, median: W1 / W2 / W3
     - 0.044 to 0.046 s / 0.435 to 0.451 s / 0.029 to 0.030 s
   * - The same computation in a fresh interpreter, median
     - 0.042 to 0.044 s / 0.432 to 0.464 s / 0.027 to 0.029 s
   * - A fresh interpreter's wall time for the same job
     - 0.60 to 0.64 s / 0.98 to 1.06 s / 0.58 to 0.62 s
   * - The template's resident size, start to end
     - 98.7 to 99.5 MB, to 155.5 to 156.5 MB
   * - Its open descriptors, unreaped children and threads at the end
     - 4, 0 and 1 in every cell
   * - Positive controls
     - the thread counter reads 1 at start, 2 with a started thread and
       2 after a BLAS call; a template that calls BLAS before forking:
       20 of 20 children killed by a signal

With 0 failures in 3 000 W1 forks, the rule of three bounds the
per-fork failure rate of a cell below 1e-3 at 95 % confidence; the bound
is per cell and per workload, the population the specification set. A
child's compute time matches the fresh interpreter's within the spread of
the cells, so the forked child computes no slower; what it saves is the
rest of the fresh wall time, about 0.55 s per miss on W1.

**The environment a generation pins is the one it was handed.** Each job
now carries the environment the asking process hands the generation
(``_Job.environment``, the projection ``_child_environment`` onto the
declared variables), and the generating process records that, not its
own ``os.environ``. The two differ: under ``LANG=C`` with no
``LC_CTYPE`` the interpreter applies PEP 538 locale coercion at its
start and adds ``LC_CTYPE=C.UTF-8`` to its own environment, so the
environment pin of every entry generated under ``LANG=C`` named a
variable the asking process never had, and every such entry was born
``Stale`` and its call raised "does not validate" (`[M]` 2026-10-16,
``ta/fork/p_lang_c.py``). It held on the fresh route before E as well.
Gate ``test_the_environment_pin_is_what_the_generation_was_handed``
(fresh and forked) generates under ``LANG=C`` and requires a ``Hit``.

**What a miss costs with the template.** `[M]` 2026-10-16, the probe
``scratch/characteristic_architecture/p3/m602/e_miss_cost.py`` (one
one-region sphere block of degree 3, twelve salted misses in one process,
each on a temporary store), run under ``-O`` in five separate processes per
route, the machine otherwise idle; the log is ``e_miss_cost.log`` beside it.

========================  ===========================  ==================
route                     median of misses 2 to 12      first miss
========================  ===========================  ==================
forked from the template  0.293 to 0.301 s              0.81 to 1.01 s
fresh interpreter         0.650 to 0.663 s              0.66 to 0.70 s
========================  ===========================  ==================

A cheap miss costs less than half of what it costs through a fresh
interpreter, because the interpreter's start and the imports are paid once
per template instead of once per miss. The first miss of a process costs
more than a fresh one, since it also starts the template and the template
builds the manifest of its own imports; every later miss repays it. All
twelve misses of every forked run were forked (``FORKS 12``) and no
template refused (``REFUSED None``).


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
  arrays (each refused when keyed);
* does not depend on a set's iteration order, which follows the hash
  seed: the key orders a set by its elements, and the generation iterates
  under ``PYTHONHASHSEED=0``, not the caller's seed (added with #592,
  2026-10-10; :ref:`verification-reference-cache-key`).

The fourth clause is enforced: each excluded argument is refused when the
key is taken. The other four are a contract, not a gate. The fifth is the
reason the key may sort a set: a function that obeys it computes the
same value whichever order its sets iterate in, so two calls whose sets
differ only in order are one call. A function that broke it would be
served the value computed under seed 0, which an in-process call under
another seed need not reproduce. The contract is
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


.. _verification-reference-cache-unseen:

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
* **Code of a source file the loader did not compile.** In a generating
  process an ``exec`` of a first-party file read by hand, or ``runpy``,
  produces code whose bytes no loader recorded; a re-read after the run
  could pin bytes that never ran, so the manifest refuses it as
  :class:`~orpheus.numerics.traced_memo.Unpinnable` (gate P7). No
  generator in ``orpheus/`` does it (`[M]` the census in
  :ref:`verification-reference-cache-pinned-source-loader`).
* :func:`~orpheus.numerics.traced_memo.trace_call` **keeps the re-read,
  and with it the window.** The in-process instrument has no loader: the
  modules it runs were imported before it was called, possibly from a
  timestamp ``.pyc``, so it pins a source file by reading it after the
  call (gate P7b), and a same-size edit within one second of the previous
  write, or an edit during the call, can be pinned beside an answer
  computed from the text before it. Its memo pin is the same read: the
  digest of the memo's three source files read after the call, since no
  loader recorded the bytes the memo ran from. Both reads go through the
  process's kept reads (:ref:`verification-reference-cache-kept-reads`),
  which serve an earlier read only while every file it looked at shows
  the same identity on disk, and keep nothing written within the last
  2 s, so they add no window of their own (`[R]`). It writes no entry and no
  source index; it is the instrument the gates measure the recorder with.
* **The modules imported before the loader.** The boot script and what
  the interpreter's start imported (``site`` runs an editable install's
  finder module, a file the installer generates, which the manifest pins
  by its bytes) were compiled before the loader existed. Their bytes are
  read when the loader is installed, the closest the process can come to
  the moment they were compiled. The boot script runs as ``__main__``,
  which CPython compiles from its source and never from a ``.pyc``, so
  its exposure is the instant between its compile and that read `[R]`.
  The finder module is imported like any other module, so a same-size
  rewrite of it within the second of its previous write leaves the
  stale ``.pyc`` hazard of ERR-115 open for it `[R]`; it is rewritten
  only by a reinstall.
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
  :math:`T_{00} = P_{ss}` identity (``greens_function_slab.py``, since
  step (e1a) under ``characteristic/origins/specular/``). A global
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

**The suite, before and after**, from the specification (the files
named ``test_trajectory_resolvent_*`` and ``test_peierls_greens_function_*``
below were deleted with their family on 2026-10-10; the tables record
P3) (`[M]` 2026-10-04, ``0e8740ef`` and the prototype, one run per file, concurrent
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
real-tree witnesses (M4.4) copied every ``*.py`` under ``orpheus/`` into a
temporary directory and made a behaviour-neutral edit there: an edit to a
construction helper the reading ran made the reading stale and left the
solve valid; an edit to the solver made both stale, the reading through
its child pin; an edit to the cylinder's rays, which a sphere reading
never ran, left both valid; a constant in the reading module's skeleton
made the reading stale. Since P1 step (e2) the witnesses are re-posed on
the characteristic reading (``test_m4_4c_…`` in
``tests/gates/numerics/test_traced_memo_clients.py``): an edit to the
wall reader ``_wall_at_face`` makes the reading stale, an edit to the polar
rule of another chart leaves it valid, and a constant in the reading
module makes it stale.

**The cadence: the specification's timing protocol** (`[M]` 2026-10-04,
``main`` at ``2bae331a``, a detached worktree on the host ``.venv``,
``python -O -m pytest --durations=0`` per file, serially; the scripts and
logs are in ``scratch/reference_architecture/p3/timing/``). Each file ran
four times, in this order over the whole list: with the memo bypassed (a
pytest plugin entering :func:`~orpheus.numerics.traced_memo.bypass` for
the session: today's arithmetic), cold (an empty store), and warm twice.
"Cold" is cold for the first file only: a later file of the cold pass
meets the entries the earlier ones wrote. The warm column is the minimum
of the two warm runs, which agreed within 2 % for every file. The machine
was not otherwise idle (a desktop session, load average about 8), so
every figure is an upper bound.

.. list-table:: Wall time per file, seconds
   :header-rows: 1
   :widths: 46 13 13 13 15

   * - File
     - Bypassed
     - Cold
     - Warm
     - Bypassed / warm
   * - ``test_peierls_greens_function_cylinder_mr.py``
     - 546.2
     - 567.2
     - 2.6
     - 210
   * - ``test_peierls_greens_function_mr.py``
     - 161.8
     - 168.9
     - 2.4
     - 67
   * - ``test_phase_c_crosscheck.py``
     - 1112.3
     - 1464.3
     - 43.2
     - 26
   * - ``test_unified_matvec_cylinder.py``
     - 848.3
     - 57.1
     - 54.1
     - 16
   * - ``test_peierls_greens_function_cylinder_mr_xverif.py``
     - 23.0
     - 26.7
     - 2.4
     - 9.6
   * - ``test_reference_body.py``
     - 19.1
     - 17.3
     - 2.8
     - 6.8
   * - ``test_trajectory_resolvent_reference.py``
     - 410.8
     - 570.8
     - 161.5
     - 2.5
   * - ``test_l1_standoff_slab_cylinder.py``
     - 1845.7
     - 1080.6
     - 1049.0
     - 1.8
   * - ``test_aba_specification.py``, ``test_crosscheck_harness.py``,
       ``test_peierls_greens_function_garcia2021.py``
     - 11.1
     - 10.5
     - 10.3
     - 1.1
   * - **The eleven files**
     - **4978.3**
     - **3963.4**
     - **1328.3**
     - **3.7**

The files that read a multi-region solve or a trajectory reading drop to
seconds once warm. What remains warm is not reference work, and it is
not meant to be memoised:

* ``test_l1_standoff_slab_cylinder.py`` spends its time in the PRODUCTION
  solves it verifies (the slab Krylov solve at 160 cells, 310 s; the
  cylinder refinement at 80 cells, 188 s). A fine-mesh production run is
  never a reference (the plan's reference ban), so these run every time;
  the memo removed the A|B|A cylinder reference those rows compare with,
  which is why the file is faster cold than bypassed.
* ``test_trajectory_resolvent_reference.py`` (deleted since) kept 91 s in the two rows
  that count solves under ``bypass()`` (a spy counts only in-process) and
  47 s in the row that recomputes the cylinder reading against a fine
  angular rule, its own independent check.
* A cold run pays one interpreter start per generation and a validation
  per hit, so a file whose references are not yet in the store
  (``test_phase_c_crosscheck.py`` first, or the cylinder solver tests) is
  slower cold than bypassed, by up to a third.

The store after the cold pass: 9 528 kB under ``.cache/references`` (``du
-sk``), 165 entries: 134 readings, 21 cylinder solves and 10 sphere solves.
The largest entry is a cylinder solve, 229 126 bytes.


.. _verification-reference-cache-limits:

Declared limits
===============

Each limit is a phase of #405 in the plan of record
(``.claude/plans/reference_cache.md``):

* **One family only.** No reference family besides the characteristic
  reference is memoised: its reading, its solve as a child entry (#592,
  ``23e7711e``) and its transport blocks as children of the solve (#405
  P3 step 2c); the other families migrate in
  phase P4, and a certification run becomes a client there. No gate yet refuses a
  traced memo on a withdrawn generator (the generating process refuses it,
  since its declared environment holds no ``ORPHEUS_*`` switch).
* **No cache on continuous integration.** The cache lives under
  ``.cache/references/`` on one machine; the continuous-integration job,
  the read-only test shards and the scheduled cold rebuild that would
  control the recorder itself are phase P5. A coarse cache key is safe
  there because a restored stale entry is a miss, never a wrong answer.
* **The suite timings are single runs.** The three files' timings
  (:ref:`verification-reference-cache-kept-reads`) are one serial run per
  tree and per mode, and the ``53551868`` runs shared the host with
  other agents' work; the protocol of the specification's section 3
  (cold, warm and bypassed, serially, three runs each) is owed.
* **Nothing prunes the memo's three stores.** An entry of the reference
  store whose pins no longer validate is regenerated under the same name,
  but an entry no key will ask for again (an argument no caller passes any
  more) stays; every source text any generation has compiled keeps its
  entry under ``.cache/bytecode/``, one per text, per interpreter and per
  ``-O`` level; and every source text any generation has indexed keeps its
  index under ``.cache/sources/``, one per text, per interpreter and per
  version of the memo's source, so every edit to the memo orphans every
  index. All three grow without bound until removed by hand: #603. `[M]`
  2026-10-15, as filed there: after the per-group block memo was renamed
  (``transport_block`` to ``group_transport``) the store held 261 entries
  under the old function id that nothing reads, and the bytecode cache
  held 366 entries after one timing run; qa's round-two review counted
  1 494 index files, 128 of them under the memo's current source. `[M]`
  2026-10-16, this checkout (``archivist_numbers.py``): 1 926 index files
  (5.1 MB), 198 under the memo's current source; 1 899 bytecode entries
  (9.7 MB); 1 261 reference entries (80.3 MB). None of the growth can
  serve a wrong value: a stale entry reads ``Stale``, a bytecode entry
  runs only under its own key, and an index is read only under the name
  of the text, interpreter and memo it was built from.
* **Forking is measured on one platform.** A miss is forked from a
  template on darwin/arm64 with Accelerate only, where the fork-safety
  protocol ran; every other platform, continuous integration's Linux with
  OpenBLAS included, generates in a fresh interpreter per miss, at the
  pre-E cost (``FORK_MEASURED``, a declared scope boundary;
  :ref:`verification-reference-cache-fork-template`).
* **A data file rewritten and restored within one run is not seen**, and
  a run that changes a file or a flag it read is refused, never cached
  (:ref:`verification-reference-cache-manifest`, "Every input at its
  first touch").
* **An** ``*.egg-info`` **install is refused.** A site-packages file is
  pinned by the version and ``RECORD`` of the distribution whose ``RECORD`` lists it,
  and an egg-info install has no ``RECORD``, so a run that touches one of
  its files is :class:`~orpheus.numerics.traced_memo.Unpinnable` and
  writes no entry. Before ``0c569df8`` such a file was pinned through the
  top-level map; the refusal is loud, never a stale hit, and `[M]`
  2026-10-15 no egg-info install is in the environment. A generator that
  needs one reopens the ruling.


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
     - Fixed: set order kept (replaced on 2026-10-10, #592, by the
       order of the elements' encodings: the kept order followed the hash
       seed); sparse class keyed and a ``dok_matrix``
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
``.claude/plans/reference_p3_spec.md``, and of what later work changed in
the memo. The 2026-10-04 rows sat on the branch
``refactor/p3-ambient-state`` at the time of writing;
``git merge-base --is-ancestor <hash> main`` outranks this column.

.. list-table::
   :header-rows: 1
   :widths: 11 55 10 24

   * - When
     - Milestone
     - Issue
     - Where
   * - 2026-10-17
     - **The 24 files re-timed after concepts 1 to 3** (the protocol of
       2026-10-14, at ``a456a9bf``). warm1 9 307 s to 3 073 s, warm2
       991 s to 487 s; system, albedos, convergence and the independent
       references take 3 to 6 s warm. 2 493 s of the remaining warm cost
       is the assembly file's 21 ``slow`` block theorems, which compute the
       system under test in process; their levers are #587 and ruling 15.
     - #405, #587
     - (this commit)
   * - 2026-10-17
     - **The PS-1982 vacuum-sphere solver becomes a client** (#405 P3
       concept 3; the user's ruling, "Memo the solver").
       ``solve_ps1982_vacuum_sphere`` carries ``@traced_memo`` alone and
       the test helper's ``lru_cache`` is removed; a memo on the helper
       was refuted (it would pin the characteristic reference), and a
       memoised eigenpair (the four cases are two eigenproblems; about
       36 s cold against 72 s) was not taken, since it restructures the
       reference. The nine rows of its gate file: 87 s cold, 2.25 s warm;
       under ``bypass()`` a 1.001 scaling of the operator moves
       :math:`k` by :math:`1.45 \times 10^{-3}`, and the store serves the
       unmutated value. A later correction refutes the review's reading
       that a docstring edit regenerated the entries: the pins strip
       docstrings, and both docstring edits left the file warm.
     - #405
     - ``da043d1c``, ``dcf631dd``, (this commit)
   * - 2026-10-17
     - **The slow cylinder point readings read through the door** (#405
       P3 concept 2; the user's ruling, "Door, with batteries in bypass").
       Eight ``slow`` rows read ``CharacteristicDerivation.evaluate(PointValue)``
       instead of ``GalerkinSystem.point_flux`` in process, through
       ``_given_by_the_door``; 13 of 13 group readings within 1 ulp; a
       battery reaches them only under ``bypass()``; C6 stays in process.
       Warm: both gate files 217 s, the 23 ``slow`` rows 25.8 s, each door
       row 0.03 to 0.13 s.
     - #405, #608, #609
     - ``94143ecd``
   * - 2026-10-17
     - **Two more oracle routes stored** (#405 P3). The partial-mirror
       block total ``specular_cylinder_total`` carries ``@traced_memo`` in
       place of its ``lru_cache``, and the lift's inline polar density
       became the stored column ``cylinder_polar_density``, one entry per
       column (a node costs less than a miss). 27 rows: cold 162 s, warm
       29 s; the polar row ``[wall-mirror]`` 14.3 s warm to 1.5 s.
     - #405
     - ``76d0409d``, ``712293a9``
   * - 2026-10-17
     - **The first test-oracle clients** (#405 P3, concept 1). A warm run
       of the slow lift and reading rows spent 993 s in two mpmath routes
       of ``_characteristic_point_mp``, ``phi_cylinder_vacuum`` and
       ``phi_cylinder``; both carry ``@traced_memo`` alone, their
       ``lru_cache`` removed (above a memo it serves the value from before
       a ``bypass()``). The other routes take milliseconds and are not
       stored. The oracle's independence is the AST census
       ``test_each_stored_reference_imports_nothing_but_the_memo``, its control
       the ``lines`` module, three first reds; a first version comparing
       fresh interpreters' module sets was blind (the memo alone loads 76
       ``orpheus`` modules; a function-body import never runs). Before it,
       the child resolution ``_resolve_memo`` (``7ca91c43``): a memo under
       a wrapping decorator generated by calling itself and was refused;
       gates M3.1b and M3.1c. 23 slow rows: warm 353 s to 127 s, cold
       1 075 s to 1 044 s, 17 entries to 34. Follow-ups: #605 (an entry
       pins the eagerly imported numerics package), #606 (seven copies of
       the ``_LAYER_SCRIPT`` import probe), #607 (a warnings filter does
       not reach the child).
     - #405; #605, #606, #607
     - ``7ca91c43``, ``333268fb``, (this commit)
   * - 2026-10-17
     - **A validation's reads are kept for the process** (#602, the hit
       cost). A warm run of the three files cost what no memo did (222 s
       against 221 s at ``3ff7866b``): 111.5 of 292.5 profiled seconds were
       validations, each re-reading every pinned source and re-finding
       every pinned distribution. One primitive in the boot module,
       ``FileIdentity``, ``Witness`` and ``Kept`` (stat before read; only
       settled reads kept, ``TRUSTED_AGE_NS`` = 2 s), and five kept reads
       built on it: ``_SOURCES`` (shadowing through witnessed candidates),
       ``_DISTRIBUTIONS``, ``_FILE_DIGESTS``, ``_LISTINGS``,
       ``_MEMO_DIGESTS``; ``_locate`` and ``_Installed._versions``
       retired. The elegance review unified the first version's two
       hand-written caches; qa's four unguarded claims became gates. Gates
       M2.3b (re-posed) to M2.3i. The three-file re-timing that is #602's
       last done-when: 221 s without the memo, 820 s cold, 163 s warm on
       the landed version (227, 725 and 155 s on the first, ``53551868``).
     - #602
     - ``8fbb1d61``
   * - 2026-10-16
     - **Each miss forked from a recording template** (#602, cut E). A
       template interpreter does the imports once and forks one generating
       child per miss; before each fork it re-validates the manifest of its
       imports (with the whole bytes of every source it compiled) and
       counts its OS threads; the child restores the recording's snapshot.
       Typed answers (``_Generated``, ``_Failed``, ``_Refused``,
       ``_Died``), one ``_Launch`` for both routes, a template replaced
       after 256 forks, forking turned off by a fresh template's refusal,
       and the platform boundary ``FORK_MEASURED``. The fork-safety
       protocol passed: 32 400 forks over 9 thread settings, 0 failures.
       Review fixes: every input (a file's bytes, a path's presence, a
       directory's names) observed at the run's first touch (ERR-116), the
       run's own output not observed, and a run that changes an input it
       observed refused; the environment pin is
       what the generation was handed (PEP 538 coercion under
       ``LANG=C``); a distribution pinned by its ``RECORD`` digest; a lock
       on the template's pipes; ``sys.argv`` emptied and standard input at
       ``/dev/null`` in the child. Gates E1 to E9 and the qa and N rows of
       ``test_traced_memo_template.py``.
     - #602
     - ``98dc5242``
   * - 2026-10-15
     - **Each source text is parsed once per content, on disk** (#602,
       round two, B). A text's index (the skeleton digest, every outermost
       def with its span and digest, every scope's span) is built from one
       parse and kept under ``.cache/sources/``, named by the bytes, the
       interpreter and the memo's source; only a generating process
       writes one, named by the bytes of the memo it ran, and the same
       digest is now a generation's memo pin (qa's finding F1, the
       elegance review's F4). Review fixes: the parse and the index one
       type (``_parsed_index``), span containment in place of a stored
       owner, one cache root and one atomic write for the three stores,
       ``_Installed`` holding ownership and versions, a file several
       ``RECORD`` files list refused, an unreadable index corrupt. Gates
       B1 to B19 and D2b; a cheap miss 0.66 to 0.85 s, the three-file cold
       pass 406 s.
     - #602
     - ``155fbdbd``
   * - 2026-10-15
     - **A site-packages file is owned by the distribution whose**
       ``RECORD`` **lists it** (#602, D). The top-level map pinned a
       namespace top's files to its first claimant
       (``sphinxcontrib/tikz.py`` to ``sphinxcontrib-applehelp``); one
       table read from every ``RECORD`` replaces it and the per-file
       ``RECORD`` search, 18 to 20 ms per manifest. A file no ``RECORD``
       lists, and every egg-info file, is refused. Gates D1 to D8.
     - #602
     - ``0c569df8``
   * - 2026-10-15
     - **A generation compiles from the bytes it pins** (#602, round
       two). A loader installed at ``sys.meta_path[0]`` by the boot script
       compiles every first-party module from the bytes it reads, keeps
       its bytecode in a content-keyed cache apart from ``__pycache__``,
       and records the bytes, which the manifest pins; code of a file it
       did not compile is refused. It closes ERR-115, a regenerated entry
       that held the old answer and validated (a timestamp ``.pyc``; an
       edit during a generation). A docstring-only def no longer makes
       ``_def_pin`` raise. The review fixes of ``8e02236f``: the skeleton
       from a shallow copy, one version lookup. Gates P1 to P9, P4b, P7b,
       A1, A2; a cheap miss 1.17 to 1.35 s.
     - #602
     - ``e666e58d``
   * - 2026-10-15
     - **A miss costs half as much** (#602): each touched file is parsed
       once, a def is digested only when a pin asks, a distribution's
       version is read once per manifest and per validation, and the
       child ends with ``os._exit``. A cheap miss went from 2.65 to 2.88 s
       to 1.14 to 1.32 s, and the three-file cold pass from 1 303 s to
       538 s; the pins are unchanged (section "What a miss costs").
     - #602
     - ``8e02236f``
   * - 2026-10-15
     - **The characteristic transport block becomes a client**,
       ``group_transport``, one group's block, a child entry of the
       solve; the three geometry types it keys on become content values.
       The miss cost it exposed (about 2.3 s of recorder and manifest per
       miss) is #602; no size threshold (the user's ruling).
     - #405, #602
     - ``9178ad00``
   * - 2026-10-10
     - **The trajectory-resolvent clients deleted.** P1 step (e2) of the
       characteristic reference campaign deleted the two multi-region
       solvers and the trajectory-resolvent reading with their family;
       the characteristic reading and its solve are the memo's clients,
       and the client gates (M4.1b to M4.8) were re-posed on them.
     - #405
     - ``5aa8cb88``
   * - 2026-10-10
     - **A set is keyed by its elements' order, not its iteration order**,
       and the generator contract gains the clause that a memoised
       function does not depend on a set's iteration order (qa's #592
       review: one solve keyed two ways across hash seeds). The
       characteristic derivation's solve becomes a memo of its own, a
       child entry of every reading. Two gates (S1 and its witness).
     - #592
     - ``45405f86``, ``23e7711e``
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
