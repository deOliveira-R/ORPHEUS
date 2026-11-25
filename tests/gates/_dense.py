r"""The declaration a test makes for a dense form it builds (N01, spec 6.2): the kernel of a Gram, named by its dimension.

Not a test module (no ``test_`` prefix). ``declared_kernel(matrix, dim)`` takes the dimension from its caller and holds no
tolerance and no cutoff: it is not a second home of ``FRAME_GRAM_RANK_RTOL`` (``orpheus/numerics/frame.py``), which is
the one rank rule, and a test that starts from a frame takes the frame's own declaration, ``frame_gram_kernel``, instead. It
spells the content key of a test-owned declaration once, ``("test-gram-kernel", content_digest(matrix))``, so one form
declared through the frame and through this helper has two digests: the key records provenance.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.declaration import Coord, Declaration, Spanned
from orpheus.values.content import content_digest


def declared_kernel(matrix: NDArray, dim: int) -> Declaration:
    r"""``Coord(n, ())`` for a full-rank form (``dim == 0``), else the ``Spanned`` of dimension ``dim`` keyed on ``matrix``."""
    n = int(np.shape(matrix)[0])
    if dim == 0:
        return Coord(n, ())
    frozen = np.array(matrix, dtype=float)
    frozen.setflags(write=False)
    return Spanned(n, ("test-gram-kernel", content_digest(frozen)), dim)
