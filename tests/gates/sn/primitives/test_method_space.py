r"""Tests for :class:`~orpheus.sn.mesh.method_space.SNMethodSpace`.

The method space is the pair (quadrature, trace ladder) an SN realizer reads. It names no face and holds no mesh: the face a
realization acts on is the ``FaceSite`` the realizer receives beside the space, so one realization carries one face value (two
independently supplied encodings of one orientation are ERR-041's own defect class). WF2.4 retired the ``face`` and ``mesh`` fields,
the ``minimal`` and ``for_face`` factories, the optional trace and ``inflow_indices_for_face``.

These tests pin:

* the field set is exactly ``(quadrature, trace)``, so a face or a mesh cannot ride on the space (D.14: the old rows asserted the
  retired fields' values; the successor asserts their absence, which a re-added field reddens);
* the trace is required (a trace-less space is unspellable, which replaces the old "refused loudly at realize" rows);
* the space holds the trace it was given, from which a face's half-spaces are read.
"""
from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from orpheus.numerics.face_layout import FaceLayout
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.spaces.angular_trace_space import AngularTraceSpace
from orpheus.sn.mesh.method_space import SNMethodSpace


pytestmark = pytest.mark.l0


def _trace(quad, faces=("xmin", "xmax")) -> AngularTraceSpace:
    layout = FaceLayout.from_named_shapes([(f, (quad.n_ordinates, 1)) for f in faces])
    return AngularTraceSpace.from_quadrature_and_layout(quad, layout)


def test_the_field_set_is_quadrature_and_trace():
    """The space names no face and holds no mesh, no index set and no factory for either."""
    assert [f.name for f in dataclasses.fields(SNMethodSpace)] == ["quadrature", "trace"]
    for retired in ("minimal", "for_face", "inflow_indices_for_face", "face", "mesh"):
        assert not hasattr(SNMethodSpace, retired), retired


def test_the_trace_is_required():
    """A quadrature-only space is unspellable: the realizers read every half-space off the ladder."""
    with pytest.raises(TypeError):
        SNMethodSpace(quadrature=Quadrature.gauss_legendre(4))  # type: ignore[call-arg]


def test_the_space_holds_the_trace_its_half_spaces_are_read_from():
    """Positive control: the held trace is the one given, and a face's inflow is the written-out incoming half-space."""
    quad = Quadrature.lebedev(17)
    trace = _trace(quad, ("xmin", "xmax", "ymin", "ymax"))
    space = SNMethodSpace(quadrature=quad, trace=trace)
    assert space.quadrature is quad and space.trace is trace
    np.testing.assert_array_equal(space.trace.inflow_indices_for_face("xmin"), np.flatnonzero(quad.mu_x > 1e-12))
    np.testing.assert_array_equal(space.trace.inflow_indices_for_face("ymax"), np.flatnonzero(quad.mu_y < -1e-12))
