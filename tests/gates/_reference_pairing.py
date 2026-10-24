r"""The ONE translation from a reference's wall to a production boundary condition, on the test side (wave 2, S5; plan 2.9 and 6.5).

Not a test module. A reference specification states its exterior as plain data, ``exterior: Mapping[int, WallAmplitudes]``:
for each tag a face carries, the cone coordinates :math:`(s, d)` of the wall's linear part (``specular``, ``diffuse``), both reals
in :math:`[0, 1]` with :math:`s + d \le 1`. A production method's problem is posed with ``boundary_conditions={tag: condition}``.
The pairing tests (a reference solution against a production solution of the SAME problem) need the second built from the first,
and that loop lives here and nowhere else:

* production does not know the reference: ``WallCondition.from_amplitudes(specular, diffuse)`` takes two floats and imports
  neither ``specification`` nor ``derivations`` (an AST gate of the conditions package);
* the reference does not know production: ``orpheus.derivations`` imports nothing from ``orpheus.transport`` or from
  ``orpheus.geometry.boundary`` (``test_references_import_no_production.py``);
* the translation therefore belongs to neither: it is this test-side loop, so the two readers of one ``WallAmplitudes`` stay
  independent derivations from one input (``instrument-doctrine`` X4), and a pairing test compares the two SOLUTIONS, never the
  translation with its own image (a test that asserted ``from_amplitudes(w)`` "realizes the same albedo as" ``Wall(w)`` would be the
  tautology, and is not written).

The key set is preserved exactly: a tag in, the same tag out, because the production exact cover refuses a condition for a tag no
face carries (``UnusedConditionError``) and a face with no condition.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING

from orpheus.transport.boundary_conditions import BoundaryConditions, WallCondition

if TYPE_CHECKING:
    from orpheus.specification.exterior import WallAmplitudes

__all__ = ["conditions_of"]


def conditions_of(exterior: Mapping[int, "WallAmplitudes"]) -> BoundaryConditions:
    """The production ``BoundaryConditions`` of a reference's ``exterior``: ``{tag: WallCondition.from_amplitudes(s, d)}``."""
    return BoundaryConditions({tag: WallCondition.from_amplitudes(w.specular, w.diffuse) for tag, w in exterior.items()})
