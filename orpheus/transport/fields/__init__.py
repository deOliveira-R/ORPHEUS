r"""L2 typed transport fields — every class inherits from
:class:`~orpheus.numerics.field.Field`.

The Depth-B unification (plan §3.3): every typed transport field is
``(values, space) + algebra`` (the L1 ABC) plus domain-specific
fields (``mesh``, ``boundary``, …). The dunder algebra is inherited
verbatim; concrete subclasses add construction factories and
diagnostics tailored to their physical role.

Migration sequence (Depth B plan §6):

* D-D — :class:`ScalarFlux` (simplest case; SHIPPED)
* D-E — :class:`HarmonicMomentField` (cleanest gap; pending)
* D-F — :class:`IsotropicSource`, :class:`PerOrdinateSource` (pending)
* D-G — :class:`BoundaryFlux` + :class:`BoundaryFaceFlux` (pending)
* D-H — :class:`AngularFlux` (most complex; pending)

Until each step lands, the corresponding ``orpheus/sn/*.py`` module
remains the canonical home for that field. A one-cycle re-export
shim is installed at the SN location AFTER each step migrates, so
in-flight callers have a window to update imports (per
``feedback_aggressive_retirement``).
"""

from __future__ import annotations

from orpheus.transport.fields.harmonic_moment_field import HarmonicMomentField
from orpheus.transport.fields.scalar_flux import ScalarFlux

__all__ = ["HarmonicMomentField", "ScalarFlux"]
