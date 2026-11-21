"""The rows of unit ``dead``, sub-unit dead-E (definitions orphaned by the W0 retirements), of the ontology W0 specification.

``Retired``: 6 rows. Each candidate was retired because a census over ``orpheus/``, ``tests/``, ``tools/``, ``examples/`` and
``student_resources/`` (an ``ast`` walk) found no production reader: ``Mixture.in_scattering_xs`` (its only consumer,
``out_scattering_xs``, retired in dead-B); ``_MMHalfGrid.downstream`` (test callers only); the ``at_order`` verb of
``GalerkinFrame`` and of the truncated bases (its only production caller, ``truncated_head``, retired in dead-D, and the
frame's own verb read the basis's). ``TransferKernel.at_order`` and ``TransferMaterialField.at_order`` are homonyms, a
different verb with production readers, and are not touched. ``KEPT``: none.
"""

from __future__ import annotations

from tests.gates._retired import Kept, Moved, Retired, Row

ROWS: tuple[Row, ...] = (
    Retired('orpheus/data/macro_xs/mixture.py', 'Mixture.in_scattering_xs'),
    Retired('orpheus/sn/angular/closure.py', '_MMHalfGrid.downstream'),
    Retired('orpheus/numerics/frame.py', 'GalerkinFrame.at_order'),
    Retired('orpheus/numerics/basis/base.py', 'TruncatedBasis.at_order'),
    Retired('orpheus/numerics/basis/legendre_basis.py', 'LegendreBasis.at_order'),
    Retired('orpheus/numerics/basis/spherical_harmonic_basis.py', 'SphericalHarmonicBasis.at_order'),
)

KEPT: tuple[Kept, ...] = ()
