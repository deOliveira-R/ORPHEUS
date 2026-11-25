r"""G03g: the form each axis carries has the bytes of the weights the axis carried before the change (N03, spec 7.4).

The record ``axis_forms/axis_weights.npz`` was taken on the tree of ``ff305e0f``, before any production edit of the unit, from the
population of ``_axis_form_population.py`` (the ``Axis.of_quadrature`` axis of each of 61 directional rules, the ``FromBasis`` axes of
the spherical-harmonic, Legendre, interval-Legendre and tensor-product families, one family with no Gram, and one synthetic measure
with a ``-0.0`` weight). The tip stores ``None`` for an all-ones axis, so the record states it as ones of the axis shape tagged
``counting``, and a family with no Gram as an empty array tagged ``noform``. This gate reads the form of the same axes on the built
tree: ``DiagonalMetric.weights`` bytes against the record, ``NoForm`` where the record says ``noform``, and a ``-0.0`` read as ``+0.0``.

The record lives beside, not inside, ``snapshots/``: the harness' loader reads every ``*.npz`` of ``snapshots/`` and compares each key.

Host-keyed: the weights of the quadrature rules pass through LAPACK, so the gate compares stored bytes on the stamped host and FAILS on a
foreign stamp with both stamps in the message (a skip would be a silent green). It is not in the CI step. Every assertion is here.

First red at the base tip: ``Axis.form`` does not exist (``AttributeError``), so the gate reads red where the spec said "none"; it is a
bit-identity witness whose behaviour is read only after the unit. Mutations: the form built from ``float32`` weights reddens through
``level_symmetric(4)``; dropping the ``-0.0`` canonicalisation of ``DiagonalMetric`` reddens through the direct leg of the second test only (the
measure's own constructor already stores ``+0.0``, so the synthetic measure's key reads green under that mutation: measured on the prototype).
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.gates.cross_method.w1 import _axis_form_population as pop
from tests.gates.cross_method.w1 import _snapshot as snap

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _require_stamped_host() -> None:
    stamp = tuple(pop.STAMP_FILE.read_text().splitlines()[0].split("\t"))
    live = snap.host_stamp()
    _require(
        stamp == live,
        f"the axis-form record cannot be read on this host: it was taken on {stamp} and this host is {live}. "
        f"Re-take it on this host from the last green commit before comparing.",
    )


def test_form_weights_bytes() -> None:
    from orpheus.numerics import metric as M

    _require_stamped_host()
    record = np.load(pop.RECORD_FILE)
    axes = pop.population()
    keys = sorted(k for k, _ in axes)
    _require(len(axes) == 71, f"the population holds {len(axes)} axes, expected 71")
    _require(sorted(k[2:] for k in record.files if k.startswith("w.")) == keys, "the record's keys are not the population's")
    bad: list[str] = []
    counts = {"weights": 0, "counting": 0, "noform": 0}
    for key, axis in axes:
        tag = str(record[f"tag.{key}"])
        stored = record[f"w.{key}"]
        counts[tag] += 1
        form = axis.form  # AttributeError on the base tip
        if tag == "noform":
            if not isinstance(form, M.NoForm):
                bad.append(f"[{key}] the record says no Gram; the form is {type(form).__name__}")
            continue
        if not isinstance(form, M.DiagonalMetric):
            bad.append(f"[{key}] the form is {type(form).__name__}, not a DiagonalMetric")
            continue
        live = np.ascontiguousarray(np.asarray(form.weights))
        if live.dtype != stored.dtype or live.shape != stored.shape or live.tobytes() != stored.tobytes():
            bad.append(f"[{key}] bytes differ from the record (tag {tag}): dtype {live.dtype} against {stored.dtype}, shape {live.shape} against {stored.shape}")
    _require(not bad, f"{len(bad)} of {len(axes)} axes differ from the record: {bad[:6]}")
    _require(counts == {"weights": 69, "counting": 1, "noform": 1}, f"the record's tags are {counts}")


def test_negative_zero_reads_as_positive_zero_in_the_form() -> None:
    r"""A ``-0.0`` weight is one value with one byte pattern: the form of the synthetic measure and a bare ``DiagonalMetric`` carry ``+0.0``."""
    from orpheus.numerics import metric as M
    from orpheus.values.content import content_digest

    key, axis = [(k, a) for k, a in pop.population() if k == "measure/negative-zero weight"][0]
    w = np.asarray(axis.form.weights)
    _require(not bool(np.signbit(w).any()), f"[{key}] the axis form keeps a negative zero: {w.tobytes().hex()}")
    d = M.DiagonalMetric(np.array([1.0, -0.0, 2.0]))
    _require(not bool(np.signbit(np.asarray(d.weights)).any()), "DiagonalMetric keeps a negative zero")
    _require(content_digest(d) == content_digest(M.DiagonalMetric(np.array([1.0, 0.0, 2.0]))), "-0.0 and +0.0 are two identities of one diagonal form")
