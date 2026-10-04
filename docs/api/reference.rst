Reference solutions (``reference``)
===================================

.. automodule:: orpheus.reference

The :mod:`orpheus.reference` package holds what a reference solution says
about the exact answer of its own equation, and how that saying is held to
account. It has six modules:

* :mod:`orpheus.reference.reading`: a reference's reading of one
  observable, :data:`~orpheus.reference.reading.ReferenceReading`, is an
  :class:`~orpheus.numerics.enclosure.Enclosure` (a value and a guaranteed
  distance to the exact one, which the reference derives), a
  :class:`~orpheus.reference.reading.Printed` value (a cited author's
  decimal text and the place it is printed, never recomputed), or an
  :class:`~orpheus.reference.reading.Uncertified` value (the reference's
  value where its family cannot yet derive a bound: no guarantee). A
  production method's reading is a different sum, so the two cannot be
  passed for each other;
* :mod:`orpheus.reference.withdrawal`: the
  :class:`~orpheus.reference.withdrawal.Withdrawal` ``(reason, issue)``, a
  maintainer ruling that a reference must not be believed until the issue
  that records it closes. It is the one class read by the test harness (the
  ``@pytest.mark.withdrawn`` marker), by the run-time lock on withdrawn
  generators (:mod:`orpheus.derivations.common.withdrawal`) and by the
  standing of the two values below;
* :mod:`orpheus.reference.published`: a
  :class:`~orpheus.reference.published.PublishedSolution`, the values a
  cited work prints for a specification, per observable, each with its own
  citation, and its :data:`~orpheus.reference.published.Standing`,
  :class:`~orpheus.reference.published.Current` or a ``Withdrawal``. It
  reads only what it prints; asked for anything else it raises
  :class:`~orpheus.reference.published.NotPrinted`;
* :mod:`orpheus.reference.certificate`: the
  :class:`~orpheus.reference.certificate.ReferenceCertificate`, a
  reference's :class:`~orpheus.reference.certificate.Claim` per observable
  (the target its bound must meet, and how its enclosure is established:
  :class:`~orpheus.reference.certificate.Exact` or
  :class:`~orpheus.reference.certificate.DerivedBound`, and nothing else),
  the evidence that can refute the claims (a published solution as an
  anchor, :class:`~orpheus.reference.certificate.Corroboration`, and a
  :class:`~orpheus.reference.certificate.Refinement` sequence, falsifying
  evidence only), and its standing. Its
  :attr:`~orpheus.reference.certificate.ReferenceCertificate.state` is
  derived on every read, never stored:
  :class:`~orpheus.reference.certificate.Valid`,
  :class:`~orpheus.reference.certificate.Invalid` with every failed check's
  reason, or the ``Withdrawal`` itself;
* :mod:`orpheus.reference.solution`: the
  :class:`~orpheus.reference.solution.ReferenceSolution`, a reference
  method's answer to a specification: its
  :class:`~orpheus.reference.solution.Derivation` evaluates any
  admissible observable on demand, the reference's natural extension with its
  derived bound where the family has one, and its certificate, when the
  family has one, claims and corroborates the observables it declares (a
  claim that disagrees with the derivation cannot be built, nor one the
  derivation leaves uncertified). Where a family cannot derive a bound the
  reading is :class:`~orpheus.reference.reading.Uncertified`. Every value but the
  reference solution itself is a content value. The reference solution is
  not itself cached: since phase P3 a derivation's ``evaluate`` may be a
  traced memo keyed on the derivation's content and the observable, whose
  entry is valid for as long as the code that ran to produce it is
  unchanged (the trajectory resolvent's is; the exact infinite medium's,
  cheaper than one interpreter start, is not), and the certificate is
  re-checked against that evaluation where the solution is built
  (:ref:`verification-reference-cache`);
* :mod:`orpheus.reference.verification`: a production answer's reading put
  on trial against a ``Valid`` reference. A
  :class:`~orpheus.reference.verification.VerificationCertificate` reads both
  sides itself and decides, exactly, the agreement
  :math:`|m - v_{\rm ref}| + b_{\rm ref} \le \tau` with the floor
  :math:`b_{\rm ref} \le \tau/10`; an
  :class:`~orpheus.reference.verification.OrderVerification` holds
  production's observed order over a refinement as intervals (each error known
  within the reference bound plus the answer's established algebraic error)
  against a declared order and band. Both refuse an uncertified reading;
  :class:`~orpheus.reference.verification.UncertifiedComparison` is the
  explicit, weaker claim a test makes against one (production within the
  tolerance of the reference's value; no verification) on exactly the
  readings the verbs cannot anchor on, and it refuses a reading a ``Valid``
  reference encloses. All three are returned, never stored.

The observables a reading is keyed on are
:mod:`orpheus.numerics.observable`; whether an observable fits a problem is
decided once, by
:func:`~orpheus.specification.specification.admit_observable`.

The package is input-tier: it imports :mod:`orpheus.specification`,
:mod:`orpheus.numerics`, :mod:`orpheus.data` and :mod:`orpheus.geometry`
and nothing above them, and :mod:`orpheus.derivations` imports it, never
the converse (:ref:`architecture-layering`). Every value in it is a
:class:`~orpheus.numerics.content.ContentIdentity`, admitted at
construction, so two spellings of one value are one value and can key a
cache (the reference solution excepted: it is rebuilt where it is
asked for, and what phase P3 caches is its derivation's evaluation). The
theory, the design's reasons and the error ontology it rests on are
:ref:`verification-reference-architecture`; the rulings are the plan of record
``.claude/plans/reference_cache.md`` (issue #405, phase P2), with the
specification ``.claude/plans/reference_p2_spec.md``; the gates are under
``tests/gates/reference/``.

Readings
--------

.. automodule:: orpheus.reference.reading
   :members:

Withdrawals
-----------

.. automodule:: orpheus.reference.withdrawal
   :members:

Published solutions
-------------------

.. automodule:: orpheus.reference.published
   :members:

The reference certificate
-------------------------

.. automodule:: orpheus.reference.certificate
   :members:

The reference solution
----------------------

.. automodule:: orpheus.reference.solution
   :members:

Verification
------------

.. automodule:: orpheus.reference.verification
   :members:
