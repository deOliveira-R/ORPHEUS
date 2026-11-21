"""The rows of unit ``defects`` of the ontology W0 specification.

``Retired``: 1 row (``MomentField._check_partner``, deleted; the base ``Field`` gate carries the check). ``KEPT``: none.
"""

from __future__ import annotations

from tests.gates._retired import Kept, Moved, Retired, Row

ROWS: tuple[Row, ...] = (
    Retired('orpheus/transport/fields/_bases.py', 'MomentField._check_partner'),
)

KEPT: tuple[Kept, ...] = ()
