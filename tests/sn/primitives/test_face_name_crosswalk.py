r"""L0 pins for the ``FaceLabel.face_name`` crosswalk (C4, issue #220).

``FaceLabel.face_name`` is the SINGLE-SOURCED rendering of the
structural face identity ``(axis_index, endpoint)`` into the
``"{axis}{min|max}"`` string world that ``FaceLayout``, the trace
space, ``SNMesh.bc``, and the sweep schedule all key on. Before C4
the crosswalk was implicit: ``boundary_face_layout`` hand-listed the
names per geometry, and the curvilinear ``"outer" → "xmax"``
translation was duplicated at two sites. These pins fix the rendering
exhaustively so a drift in the crosswalk fails HERE, at L0, not as a
mis-keyed boundary operator three layers up.

Verification design:
``.claude/agent-memory/test-architect/c4_snmesh_bc_dict_verification.md``
(G1.1 exhaustive table, G1.2 fail-loud negative, G2.14 d=3 synthetic
admission — the crosswalk is a pure function, so d=3 needs no Mesh3D).

The expectation table below is hand-transcribed (mirror-not-import,
vv L11): it must NOT be generated from ``AXIS_NAMES`` or
``_ENDPOINT_SUFFIX``, otherwise the test is tautological against the
production derivation it verifies.
"""
from __future__ import annotations

import pytest

from orpheus.sn.axis import FaceLabel

pytestmark = [pytest.mark.foundation]


# Hand-transcribed (axis_index, endpoint) → face-name table covering
# every canonical endpoint on every named axis, d ∈ {1, 2, 3}.
# "outer" renders as the max face of its axis: a solid radial axis's
# outer surface IS its high face (the historical curvilinear
# convention — sphere/cylinder boundary operators key on "xmax").
_EXPECTED_FACE_NAME = {
    (0, "min"): "xmin",
    (0, "max"): "xmax",
    (0, "outer"): "xmax",
    (1, "min"): "ymin",
    (1, "max"): "ymax",
    (1, "outer"): "ymax",
    (2, "min"): "zmin",
    (2, "max"): "zmax",
    (2, "outer"): "zmax",
}


class TestFaceNameCrosswalk:
    def test_exhaustive_axis_endpoint_table(self) -> None:
        """Every canonical (axis, endpoint) pair renders per the table."""
        for (axis_index, endpoint), expected in _EXPECTED_FACE_NAME.items():
            rendered = FaceLabel(axis_index, endpoint).face_name
            if rendered != expected:
                pytest.fail(
                    f"FaceLabel({axis_index}, {endpoint!r}).face_name = "
                    f"{rendered!r}; expected {expected!r}"
                )

    def test_d3_renders_z_faces(self) -> None:
        """The d=3 admission pin: a third axis names z-faces, today.

        No Mesh3D exists (C5); the crosswalk being a pure function on
        ``FaceLabel`` is exactly what makes the 3-axis rendering
        verifiable NOW, so ``boundary_face_layout``'s derived loop is
        z-correct by construction the day a 3-axis mesh constructs.
        """
        rendered = [
            FaceLabel(a, ep).face_name
            for a in range(3)
            for ep in ("min", "max")
        ]
        if rendered != ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"]:
            pytest.fail(f"d=3 face-name inventory wrong: {rendered}")

    def test_noncanonical_endpoint_fails_loud(self) -> None:
        """An overridden endpoint label has NO face name (ValueError).

        ``AxisMesh.label_low`` / ``label_high`` are user-overridable
        (e.g. ``"left"`` / ``"right"``); a renamed endpoint must not
        silently desynchronize from the ``"{axis}{min|max}"`` world.
        """
        with pytest.raises(ValueError, match="left"):
            FaceLabel(0, "left").face_name

    def test_unnamed_axis_fails_loud(self) -> None:
        """An axis beyond the named inventory raises (IndexError)."""
        with pytest.raises(IndexError):
            FaceLabel(len("xyz"), "min").face_name
