r"""Transport-layer discrete frames — carrier-typed specializations of the
numerics :class:`~orpheus.numerics.frame.FrameBase` hierarchy.

The generic frames in :mod:`orpheus.numerics.frame` are **carrier-agnostic**:
their analysis/reconstruction faces consume and emit bare ``np.ndarray`` (the
same generic mechanism serves the angular spherical-harmonic projection AND the
indicator-basis spatial homogenisation). The casting to the typed transport
:class:`~orpheus.numerics.field.Field` carriers cannot live in numerics — the
deepest primitive the carriers share, ``Field``, is in numerics, but the
*castability* (the carrier-minted spaces the leaves are elements of) lives in
the transport :class:`~orpheus.transport.fields._bases.BulkField` base. So the typed
seam lives HERE, in transport, one layer above numerics: the
:class:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame` factory MINTS
bound, carrier-typed operator faces
(:class:`~orpheus.transport.frames.harmonic_frame.HarmonicAnalysisOperator` /
:class:`~orpheus.transport.frames.harmonic_frame.HarmonicReconstructionOperator`
— members of the numerics analysis/reconstruction ROLE family) over the
generic faces (F-1, ``frame_square_recarve.md``).
"""

from __future__ import annotations

from orpheus.transport.frames.harmonic_frame import (
    HarmonicAnalysisOperator,
    HarmonicFrame,
    HarmonicReconstructionOperator,
)

__all__ = [
    "HarmonicAnalysisOperator",
    "HarmonicFrame",
    "HarmonicReconstructionOperator",
]
