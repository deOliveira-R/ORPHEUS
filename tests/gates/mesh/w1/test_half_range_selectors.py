r"""G02a's nudged-fixture row: every selector reads the half-range owner (N02, spec 7.3).

The fixture is the shipped ``product(2, 4)`` with its 4 exactly-zero ``mu_x`` replaced by ``9.5e-16`` (a ``dataclasses.replace`` of the
measure). For ``xmin`` and ``xmax`` the four readers of "which ordinates leave the face" must give one index set: ``face_outflow_ordinates``,
``AngularTraceSpace.outflow_indices_for_face`` and ``inflow_indices_for_face``, the x label of ``Quadrature.octants`` and
``half_range(quad, 0, sign=+-1)``. The file lives under ``tests/gates/mesh/w1`` because it calls ``face_outflow_ordinates``, and the
test-layer gate refuses an import of ``orpheus.mesh`` from ``tests/gates/numerics``.

Exact on any host: integer index sets of a fixed fixture, no LAPACK or libm result compared.

First red today: the live readers disagree with each other, before the owner is looked at: the trace reads ``[0, 1, 3, 4, 5, 7]`` (its
``4 eps`` band is below ``9.5e-16``), ``face_outflow_ordinates`` reads ``[0, 4]`` (a ``1e-15`` cutoff is above it), and the octant x
label is ``0`` on the 4 nudged ordinates. The three rebinding mutants (each consumer's ``half_range`` name replaced by a ``1e-15``
tolerance classifier) must each redden the agreement row; they can be installed only once the consumers bind the name, so they are red
(``AttributeError``: the module has no ``half_range``) until the unit is built, and the same ``4 eps`` copy does not redden it at
``9.5e-16``.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from orpheus.mesh import axis as mesh_axis
from orpheus.mesh.axis import face_outflow_ordinates
from orpheus.numerics.coproduct_layout import CoproductLayout
from orpheus.numerics.quadrature import directional
from orpheus.numerics.spaces import AngularTraceSpace, angular_trace_space
from tests.gates import _half_range_helpers as H

pytestmark = pytest.mark.foundation

_FACES = (("xmin", "min", -1), ("xmax", "max", +1))
_NG = 2


def _owner():
    from orpheus.numerics.tangential import half_range

    return half_range


def disagreements(quad) -> list[str]:
    r"""What the readers of the half-range say differently on ``quad``, face by face; an empty list is one answer. Owner not consulted."""
    n = int(quad.n_ordinates)
    layout = CoproductLayout.from_named_shapes([(f, (n, _NG)) for f, _, _ in _FACES])
    trace = AngularTraceSpace.from_quadrature_and_layout(quad, layout)
    x_label = H.octant_labels(quad)
    x_label = x_label if x_label.ndim == 1 else x_label[:, 0]
    bad: list[str] = []
    for face, endpoint, sign in _FACES:
        label = SimpleNamespace(axis_index=0, endpoint=endpoint)
        out = {
            "face_outflow_ordinates": tuple(face_outflow_ordinates((), label, quad).tolist()),
            "trace outflow": tuple(trace.outflow_indices_for_face(face).tolist()),
            "octant label": tuple(np.flatnonzero(x_label == sign).tolist()),
        }
        inn = {
            "trace inflow": tuple(trace.inflow_indices_for_face(face).tolist()),
            "octant label": tuple(np.flatnonzero(x_label == -sign).tolist()),
        }
        if len(set(out.values())) != 1:
            bad.append(f"{face} outflow: {out}")
        if len(set(inn.values())) != 1:
            bad.append(f"{face} inflow: {inn}")
    return bad


def test_every_selector_reads_the_owner_on_a_nudged_fixture() -> None:
    quad = H.nudged_product(9.5e-16)
    mu_x = np.asarray(quad.mean_axis_cosine(0))
    assert int(np.sum((mu_x > 0) & (mu_x < 1e-10))) == 4, "the fixture is not nudged (activation)"

    bad = disagreements(quad)
    assert not bad, f"the readers disagree among themselves: {bad}"

    half_range = _owner()
    n = int(quad.n_ordinates)
    layout = CoproductLayout.from_named_shapes([(f, (n, _NG)) for f, _, _ in _FACES])
    trace = AngularTraceSpace.from_quadrature_and_layout(quad, layout)
    for face, endpoint, sign in _FACES:
        half = half_range(quad, 0, sign=sign)
        assert tuple(trace.outflow_indices_for_face(face).tolist()) == tuple(half.outflow.indices)
        assert tuple(trace.inflow_indices_for_face(face).tolist()) == tuple(half.inflow.indices)
        label = SimpleNamespace(axis_index=0, endpoint=endpoint)
        assert tuple(face_outflow_ordinates((), label, quad).tolist()) == tuple(half.outflow.indices)
    xmax = half_range(quad, 0, sign=+1)
    assert (len(xmax.outflow.indices), len(xmax.inflow.indices)) == (6, 2), (xmax.outflow.indices, xmax.inflow.indices)


def test_the_agreement_instrument_reads_no_disagreement_on_the_exact_zero_fixture() -> None:
    r"""Control that runs on any tree: on the exact-zero shipped fixture the readers agree (no disagreement is read for nothing)."""
    from orpheus.numerics.quadrature import Quadrature

    assert disagreements(Quadrature.product(n_mu=2, n_phi=4)) == []


@pytest.mark.parametrize("mod", [mesh_axis, angular_trace_space, directional], ids=["mesh.axis", "angular_trace_space", "directional"])
def test_a_banded_owner_in_one_consumer_reddens_the_agreement_row(mod, monkeypatch: pytest.MonkeyPatch) -> None:
    r"""Mutation: the ``half_range`` name that the module binds replaced by the ``1e-15`` tolerance classifier; the row must read a disagreement.

    Precondition: the honest fixture reads none. The precondition and the patch both need the module to bind ``half_range``
    (``raising=True``), so each arm is red until its consumer calls the owner.
    """
    quad = H.nudged_product(9.5e-16)
    assert disagreements(quad) == [], "the readers disagree before the mutation"
    monkeypatch.setattr(mod, "half_range", H.banded_owner(1e-15))
    after = disagreements(H.nudged_product(9.5e-16))
    assert after, f"{mod.__name__}: a 1e-15 tolerance classifier in this consumer reads the same as the owner (the mutation did not bite)"


def test_the_four_eps_copy_does_not_redden_the_agreement_row_at_the_round_off_cosine(monkeypatch: pytest.MonkeyPatch) -> None:
    r"""``4 eps`` = 8.9e-16 is below ``9.5e-16``: the retired band read this cosine as genuine too, so this row is blind to it by design
    (the subnormal stub of the value row is the arm that reddens it)."""
    monkeypatch.setattr(angular_trace_space, "half_range", H.banded_owner(4.0 * H.EPS))
    assert disagreements(H.nudged_product(9.5e-16)) == []


def test_face_outflow_ordinates_refuses_a_non_finite_cosine() -> None:
    r"""Today the selector accepts ``[nan, 0.5]`` and ``[inf, -0.5]`` and returns ``[1]`` and ``[0]``; the owner refuses them."""
    label = SimpleNamespace(axis_index=0, endpoint="max")
    for cosines in ([np.nan, 0.5], [np.inf, -0.5]):
        with pytest.raises(ValueError, match="non-finite direction cosine"):
            face_outflow_ordinates((), label, H.CosineStub(cosines))
