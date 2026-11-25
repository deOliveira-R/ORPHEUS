r"""G01j: the bare constructor takes the matrix AS STORED (N01, spec 3.3, 7.2).

The asymmetric ``G+`` of ``frame_faces.angular_frame(Quadrature.folded_product(2, 4), 2)`` (asymmetry of order 1e-17, the shape of
``test_axis_identity.py::test_I5b`` and ``test_I5c``) is handed to ``DenseMetric(m, kernel=k)`` with the kernel the frame declares for
its Gram. The second face equals ``np.linalg.pinv(m, hermitian=True, rcond=1e-12)`` of the matrix AS GIVEN, byte for byte: a
computation against a computation (the asymmetry is the host's round-off where it has one, and one ulp put in by hand where it
has none). A bare path that symmetrised its input first would move these bytes (it moves 2
of the 234 bare matrices of the build's spy comparison). Activation: the matrix is asymmetric, so the symmetrisation is a change.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics import frame as FR
from orpheus.numerics import frame_faces as FF
from orpheus.numerics.metric import DenseMetric, Face
from orpheus.numerics.quadrature import Quadrature

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


@pytest.mark.rests_on("tests/gates/cross_method/w1/test_dense_inverse_bytes.py::test_bare_path_equals_the_stored_metric")
def test_bare_path_takes_the_matrix_as_stored() -> None:
    gram = np.asarray(FF.angular_frame(Quadrature.folded_product(2, 4), 2).discrete_gram, float)
    g_sym = (gram + gram.T) / 2.0
    kernel = FR.frame_gram_kernel(gram)
    m = np.linalg.pinv(g_sym, hermitian=True, rcond=1e-12)  # G+ of the frame head, as numpy builds it
    if m[0, 1] == m[1, 0]:  # a host whose pinv returns this pair exactly symmetric: one ulp makes the matrix asymmetric
        m[0, 1] = np.nextafter(m[0, 1], np.inf)
    _require(float(np.max(np.abs(m - m.T))) > 0.0, "activation: G+ must be asymmetric for the symmetrisation to be a change")
    face = np.asarray(DenseMetric(m, kernel=kernel, face=Face.METRIC).inverse_matrix)
    expected = np.linalg.pinv(m, hermitian=True, rcond=1e-12)
    _require(face.tobytes() == expected.tobytes(), f"the second face differs from pinv of the matrix as given: max|a-b|/max|b| = {float(np.max(np.abs(face - expected)) / np.max(np.abs(expected))):.3e}")
    sym = (m + m.T) / 2.0
    _require(np.linalg.pinv(sym, hermitian=True, rcond=1e-12).tobytes() != expected.tobytes(), "activation: symmetrising the input must change the bytes")
