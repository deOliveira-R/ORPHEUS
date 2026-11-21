Traced memo (``memo``)
======================

The :mod:`orpheus.memo` package sits between :mod:`orpheus.values` and ``numerics``, beside ``geometry``: it imports the
standard library, numpy and ``orpheus.values`` only, and neither ``values`` nor ``geometry`` imports it (the gate is
``tests/gates/test_layer_imports.py``; the layer is on :doc:`/architecture/layering`). It was
``orpheus.numerics.traced_memo`` and its two private modules until the move that made it a package of its own.

The traced memo — :mod:`orpheus.memo.traced_memo`, a pure function of
content-identified arguments memoised on disk under ``.cache/references/``
and keyed on what ran: on a miss the call is generated in a fresh
interpreter (or, on a platform where forking is measured, a child forked
from a template interpreter that has only imported) that records, from its first line, every function that starts,
every file it reads or probes and every directory it lists
(:mod:`orpheus.memo._traced_memo_boot`); the entry keeps that
:class:`~orpheus.memo.traced_memo.Manifest` beside an exact payload,
and a lookup validates the manifest against the checkout without importing
or running anything in its own process (the one exception: an import set
with a changed module has its import outcome measured again in a fresh
interpreter), reading again only what changed on disk since the
process last read it (:class:`~orpheus.memo._traced_memo_boot.Kept`). The key is
:func:`~orpheus.values.content.encode_exact` of the function, its
signature-bound arguments and the platform tag. The theory, with why the
identity of a generator is its recorded execution, why a miss runs in a
fresh process, each pin kind and what it guards, the payload, the bypass,
the relation to the reference certificate, the clients, the limits no
manifest can see and the review round's defects, is
:ref:`verification-reference-cache`; its second reader, which reads a
recording as pins with every refusal kept as a value
(:func:`~orpheus.memo.traced_memo.collect_pins`), is the cone store
(:ref:`verification-cone-store`).

.. automodule:: orpheus.memo.traced_memo
   :members: traced_memo, TracedMemo, bypass, cache_root, default_root,
      Hit, Absent, Stale, Corrupt, Manifest, DefPin, StatementsPin,
      ImportOutcomePin, DistributionPin, InterpreterPin, DataPin, PresencePin, ListingPin,
      EnvironmentPin, WorkingDirectoryPin, ChildPin, MemoPin, validate,
      trace_call,
      origin, function_digest, platform_tag, fork_platform,
      encode_payload, decode_payload, Unpinnable, Unencodable, TemplateRefused,
      collect_pins, CollectedPins, PinRefusal, AmbiguousOwner, UnownedSitePackage,
      UninstalledDistribution, NoDefAtLine, ChangedByRun, NoRecordedBytes,
      SourceGone, UnadmittedProcess, UnnamedModule, UnmeasurableImport

.. autodata:: orpheus.memo.traced_memo.Pin
   :no-value:

.. autodata:: orpheus.memo.traced_memo.Verdict
   :no-value:

.. automodule:: orpheus.memo._traced_memo_boot
   :members: Recording, ReadChild, note_child, spawning, main, template, FORK_MEASURED,
      FileIdentity, Witness, Kept, TRUSTED_AGE_NS, file_digest, listing_digest,
      fork_blockers

The import-time statements of a module, as the memo's ``StatementsPin``
reads them (:mod:`orpheus.memo._import_time`, private to the memo and
one of its pinned sources): each top-level and class-body statement with
its digest, the names it binds and reads and whether it is an effect; the
selection a run's entry names make of them; the closure of the entry names
over the statements. The rule, the attack rows it was built against and
what it cannot see are :ref:`verification-cone-store-statements`.

.. automodule:: orpheus.memo._import_time
   :members: Statement, string_names, code_names, span, statements, holding,
      selection, closure, observed_digests
