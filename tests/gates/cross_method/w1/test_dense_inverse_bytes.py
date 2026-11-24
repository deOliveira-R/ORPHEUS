r"""G01a: the bare constructor equals the stored bytes, and the declared dimension is the nullity (N01, spec 7.2).

The harness' G00a compares the ``inverse_of`` path of the 21 dense Grams of group (c) with the stored bytes. This file reads the
same stored bytes through the OTHER path (instrument doctrine X4: two spellings of one comparison agree tautologically): the bare
constructor ``DenseMetric(g_sym, kernel=k)`` on the symmetrised stored input Gram, whose second face (``.inverse_matrix``) must
equal the stored ``metric`` face (``G+``) byte for byte. After N01 both paths call the one body ``_exact_inverse``, so a body
mutation reddens G00a and this file alike; this file's own reach is the bare constructor's wiring of ``__post_init__`` to the body.
The second test reads the declared dimension against the nullity of the stored ``G+``, ``n - round(trace(G+ G))``.

Host-keyed: it compares live bytes with the macOS-arm64 bytes of ``snapshots/frames.npz``, so on a host whose stamp differs from
the manifest's it FAILS with both stamps in the message (a skip would be a silent green), as G00a does, and it is not listed in
the CI step. Every assertion is here: a helper's ``assert`` is removed by ``python -O``.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.metric import DenseMetric
from tests.gates.cross_method.w1 import _snapshot as snap
from tests.gates.numerics.w1 import _stored_grams as sg

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _require_stamped_host() -> None:
    stamp, _ = snap.read_manifest(snap.SNAPSHOT_DIR / "MANIFEST.txt")
    live = snap.host_stamp()
    _require(
        stamp == live,
        f"bit identity cannot be read on this host: the manifest was taken on {stamp} and this host is {live}. "
        f"Re-take the snapshots on this host at the last green commit before comparing.",
    )


@pytest.mark.parametrize("family", sg.families())
@pytest.mark.rests_on("tests/gates/cross_method/w1/test_snapshot_harness.py::test_every_snapshot_equal")
def test_bare_path_equals_the_stored_metric(family: str) -> None:
    """``DenseMetric(g_sym, kernel=k).inverse_matrix`` equals the stored ``c.<family>.metric`` bytes, on each of the 21 families.

    Input: the stored raw Gram symmetrised exactly (``(g + g.T) / 2``), the kernel the harness builder declares
    (``frame_gram_kernel`` for the 17 frame-built, the builder's own for the 4 hand-written). First red today: the ``kernel``
    keyword does not exist (``TypeError``, 21 of 21). Mutation (7.2): the recombination ``(vt.T * sinv) @ u.T`` reddens 18 of 21.
    """
    _require_stamped_host()
    e = sg.stored(family)
    kernel = sg.declared_for(e)
    m = DenseMetric(e.g_sym, kernel=kernel)
    live = np.asarray(m.inverse_matrix)
    _require(live.shape == e.metric.shape and live.dtype == e.metric.dtype, f"[{family}] shape/dtype {live.shape}/{live.dtype} against stored {e.metric.shape}/{e.metric.dtype}")
    _require(
        live.tobytes() == e.metric.tobytes(),
        f"[{family}] the bare path's second face differs from the stored G+ bytes: max|a-b|/max|b| = "
        f"{float(np.max(np.abs(live - e.metric)) / np.max(np.abs(e.metric))):.3e}",
    )


@pytest.mark.parametrize("family", sg.families())
@pytest.mark.rests_on("tests/gates/numerics/w1/test_declared_kernel_is_exact.py::test_declared_kernel_is_exact_for_the_frames")
def test_declared_dim_is_the_nullity(family: str) -> None:
    """The declared dimension equals ``n - round(trace(G+ G))`` with ``G+`` the stored ``metric`` and ``G`` the symmetrised Gram.

    The trace of the projector ``G+ G`` is the rank, read from stored bytes and not from the declaration under test. First red
    today: the module under test does not import (``frame_gram_kernel``, ``numerics.declaration``). Mutation: declare ``dim + 1``
    (red on 21 of 21), declare ``dim - 1`` (red on the 14 deficient Grams).
    """
    _require_stamped_host()
    e = sg.stored(family)
    kernel = sg.declared_for(e)
    n = e.g_sym.shape[0]
    nullity = n - int(round(float(np.trace(e.metric @ e.g_sym))))
    _require(kernel.dim == nullity, f"[{family}] declared dim {kernel.dim} against the nullity {nullity} of the stored G+ (n = {n})")


def test_population_is_the_21_families() -> None:
    """The population is not empty and not a stale subset: 21 families, 17 frame-built and 4 hand-written, 14 deficient."""
    fams = sg.families()
    _require(len(fams) == 21, f"{len(fams)} stored families, expected 21")
    _require(sum(sg.stored(f).frame_built for f in fams) == 17, "expected 17 frame-built families")
    _require(set(sg.HAND_DIMS) <= set(fams), "a hand-written family id is not stored")
    deficient = sum(
        1 for f in fams
        if sg.stored(f).g_sym.shape[0] - int(round(float(np.trace(sg.stored(f).metric @ sg.stored(f).g_sym)))) > 0
    )
    _require(deficient == 14, f"{deficient} deficient stored Grams, expected 14")
