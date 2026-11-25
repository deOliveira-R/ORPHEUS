r"""The 21 stored dense Grams of the bit-identity harness (group c), and the kernel each is declared with.

Not a test module (no ``test_`` prefix) and it asserts nothing: a helper imported by a test RETURNS data (``python -O`` strips an
``assert``). The Grams are read from ``tests/gates/cross_method/w1/snapshots/frames.npz`` (``c.<family>.input_gram``, ``.metric``
and ``.gram``), so a gate that uses them compares a computation with a computation on a stored INPUT; only a gate that reads the
stored ``metric`` bytes is host-keyed (G01a, which says so).

The declared kernel is the harness builder's own: ``frame_gram_kernel`` (LIVE, through the module attribute, so a gate's
in-process mutation of the frame's count reaches it) for the 17 frame-built Grams, and the dimension written here for the 4
hand-written ones (``syn.rank2_3x3`` has rank 2, ``syn.near_cutoff_5x5`` has one value of 5e-13 below the frame's cutoff, the
other two are full rank): an independent witness of the frame's rule on those four.
"""

from __future__ import annotations

import functools
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics import frame as FF
from orpheus.numerics.declaration import Declaration
from tests.gates._dense import declared_kernel

STORED_FRAMES = Path(__file__).resolve().parents[2] / "cross_method" / "w1" / "snapshots" / "frames.npz"

#: The four hand-written families and the nullity each is declared with (the builder's own declaration, 6.2).
HAND_DIMS: dict[str, int] = {"syn.rank2_3x3": 1, "syn.pou_3x3": 0, "syn.kappa1e4_8x8": 0, "syn.near_cutoff_5x5": 1}


@dataclass(frozen=True)
class StoredGram:
    family: str
    gram: NDArray  # the raw input Gram of the builder
    metric: NDArray  # the stored ``.matrix`` face of ``DenseMetric.inverse_of(gram)``: G+
    gram_face: NDArray  # the stored ``.inverse_matrix`` face: the symmetrised Gram

    @property
    def g_sym(self) -> NDArray:
        return (self.gram + self.gram.T) / 2.0

    @property
    def frame_built(self) -> bool:
        return self.family not in HAND_DIMS


@functools.cache
def _load() -> dict[str, StoredGram]:
    with np.load(STORED_FRAMES) as z:
        names = sorted(k[len("c."):-len(".input_gram")] for k in z.files if k.endswith(".input_gram"))
        return {f: StoredGram(f, z[f"c.{f}.input_gram"], z[f"c.{f}.metric"], z[f"c.{f}.gram"]) for f in names}


def families() -> tuple[str, ...]:
    """The 21 family ids, sorted."""
    return tuple(_load())


def stored(family: str) -> StoredGram:
    return _load()[family]


def declared_for(entry: StoredGram) -> Declaration:
    """The kernel the harness builder declares for this Gram."""
    if entry.frame_built:
        return FF.frame_gram_kernel(entry.gram)
    return declared_kernel(entry.g_sym, HAND_DIMS[entry.family])
