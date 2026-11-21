"""The rows of sub-unit ``angular-mass`` (#708, ``ontology_w0_spec.md`` 5.7.3 E6): the cross-class source dunder and its deferral arm.

Five ``Retired`` rows, each DEFINED at the base tip ``954c4a4f`` (so each reads red on the real tree until the unit lands): the arm
``AngularSourceSink.__add__`` and ``__radd__`` with ``AngularSourceSink.CONTAINS``, and ``BulkField.CONTAINS`` with ``BulkField.__add__``
(whose only job was the deferral to that arm; a census with a positive control reads 3 ``CONTAINS`` nodes today, 5.7.2 item 7).
After the unit the mixed ``ScalarSourceSink + AngularSourceSink`` sum is refused by ``Field._check_partner``, the existing same-class
gate (gate G3 of ``tests/gates/transport/fields/test_scalar_source_does_not_add_to_an_angular_one.py``). ``KEPT = ()``.

Mutation (the mechanism of 1.5, applied by the reviewers): copy ``orpheus/transport/fields/_bases.py`` to a temporary root at the same
relative path, run the walk on that unmutated root and require the five rows green there (a ``Retired`` row whose file is absent is
green, so a one-file root is legal), insert ``PLANTED_BULK_FIELD_CONTAINS`` of ``_controls.py`` after the end of the ``BulkField`` class
docstring in the copy, call the walk on that root and require the ``BulkField.CONTAINS`` row in the red set, by the per-row verdict.
"""

from __future__ import annotations

from tests.gates._retired import Kept, Retired, Row

ROWS: tuple[Row, ...] = (
    Retired('orpheus/transport/fields/angular_source_sink.py', 'AngularSourceSink.__add__'),
    Retired('orpheus/transport/fields/angular_source_sink.py', 'AngularSourceSink.__radd__'),
    Retired('orpheus/transport/fields/angular_source_sink.py', 'AngularSourceSink.CONTAINS'),
    Retired('orpheus/transport/fields/_bases.py', 'BulkField.CONTAINS'),
    Retired('orpheus/transport/fields/_bases.py', 'BulkField.__add__'),
)

KEPT: tuple[Kept, ...] = ()
