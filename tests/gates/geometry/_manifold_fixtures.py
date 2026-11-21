"""Test-side fixture of the manifold catalogue: the unit interval that left production in wave W0, unit ``dead``.

``UNIT_INTERVAL`` (row 27 of the dead census) was a module value of ``orpheus/geometry/manifold.py`` that no production path read; its
five consuming test files import it from here, and it has this one home. The shipped members (``COSINE_INTERVAL``, ``HALF_LINE``,
``REAL_LINE``) stay in production.
"""

from __future__ import annotations

from orpheus.geometry.manifold import Interval

#: The interval :math:`[0, 1]`.
UNIT_INTERVAL = Interval(0.0, 1.0)
