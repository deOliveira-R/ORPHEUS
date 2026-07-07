r"""``face_streaming_normal`` — the spatial-trace partial-current measure.

``face_streaming_normal(coeff, weight) = |coeff| * weight`` computes the
magnitude of the normal streaming flux through a codim-1 **spatial**
boundary face, weighted by the angular quadrature: with ``coeff = Ω·n_f``
it is the partial-current metric ``G_s = |Ω·n_f|·w_n`` (Lewis & Miller
§3.7), the boundary inner product under which the ``BoundaryOperator``
Hilbert adjoints are self-adjoint. It vanishes at **grazing** incidence
(``Ω ⊥ n_f``).

Scope note: this is the SPATIAL-trace measure only. The angular pole
(ψ½ / starting-direction) block does NOT route through it — its Hilbert
metric is the radial cell volume ``V_cell`` (a STATE metric), not the
angular through-flux coefficient. The two codim-1 metrics do not unify
through one kernel; see
:mod:`orpheus.numerics.spaces.starting_direction_space` and ERR-067.

Foundation gates: the kernel definition, the trace-metric reproduction
at 0 ULP, the drop-abs mutation tooth, and the Mode-11 execution
sentinel (the trace build actually routes through the kernel).

vv Mode-8: ``np.testing`` / ``pytest.fail`` / ``pytest.raises`` only (the
canonical run is ``python -O``; bare asserts in helpers would be stripped).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.face_layout import FaceLayout, face_streaming_normal
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.spaces.angular_trace_space import AngularTraceSpace
import orpheus.numerics.spaces.angular_trace_space as ats_mod

pytestmark = pytest.mark.foundation


# ── independent references (NO assertions in helpers — vv Mode-8) ───────
def _inline_partial_current(space: AngularTraceSpace, quad: Quadrature) -> np.ndarray:
    """``|Ω·n_f|·w_n`` flat build spelled inline — the DEFINITION, not the SUT.

    Bit-identical to the ``_build_trace_metric_weights`` body (``np.abs(Ω·n) *
    w`` + the per-slot broadcast); it does NOT call ``face_streaming_normal``,
    so a buggy kernel diverges from it.
    """
    layout = space.layout
    ref = np.zeros((int(layout.total_size),), dtype=float)
    w = np.asarray(quad.weights, dtype=float)
    for f_idx, face in enumerate(layout.faces):
        slot = layout.faces[face]
        face_w = np.abs(space.omega_dot_n[f_idx]) * w
        face_w_axis0 = face_w.reshape((face_w.shape[0],) + (1,) * (len(slot.shape) - 1))
        flat_face = np.broadcast_to(face_w_axis0, slot.shape).reshape(-1)
        ref[slot.offset : slot.offset + slot.flat_size] = flat_face
    return ref


def _slab(quad: Quadrature, ng: int = 2) -> AngularTraceSpace:
    layout = FaceLayout.from_named_shapes(
        [("xmin", (quad.N, ng)), ("xmax", (quad.N, ng))]
    )
    return AngularTraceSpace.from_quadrature_and_layout(quad, layout)


# ════════════════════════ kernel definition ════════════════════════════
class TestKernelDefinition:
    def test_kernel_is_abs_coeff_times_weight(self):
        """|coeff|·weight — signed coeff, positive weight → positive result.

        Pins the abs on the FIRST operand (catches drop-abs AND
        abs-on-the-wrong-operand: correct(-3, 5) = +15, both mutations
        give -15).
        """
        coeff = np.array([-3.0, 2.0, -0.5, 0.0])
        weight = np.array([5.0, 7.0, 4.0, 9.0])
        np.testing.assert_array_equal(
            face_streaming_normal(coeff, weight), np.abs(coeff) * weight
        )
        np.testing.assert_array_equal(
            face_streaming_normal(np.array([-3.0]), np.array([5.0])), [15.0]
        )

    def test_measure_vanishes_at_grazing(self):
        """The measure |Ω·n|·w vanishes as the streaming-normal coefficient
        → 0 (grazing incidence, Ω ⊥ n_f) — a tangential ordinate carries no
        partial current through the face."""
        omega_dot_n = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])  # 0.0 = grazing
        w = np.ones_like(omega_dot_n)
        result = face_streaming_normal(omega_dot_n, w)
        np.testing.assert_array_equal(result, np.abs(omega_dot_n) * w)
        np.testing.assert_array_equal(result[2], 0.0)  # the grazing ordinate

    def test_one_kernel_reproduces_trace_metric_at_0ulp(self):
        """coeff = Ω·n_f reproduces the partial-current metric |Ω·n|·w bit-exact."""
        quad = Quadrature.gauss_legendre(8)
        space = _slab(quad)
        np.testing.assert_array_equal(
            space.inner_product_weights, _inline_partial_current(space, quad)
        )


# ════════════════════════ mutation tooth (in-process monkeypatch) ══════
class TestMutationTeeth:
    def test_drop_abs_reds_trace_metric(self, monkeypatch):
        """Drop the abs → signed metric; inflow ordinates (Ω·n<0) flip sign.
        The trace gate MUST diverge from the |Ω·n|·w reference and go
        negative."""
        monkeypatch.setattr(
            ats_mod, "face_streaming_normal", lambda c, w: np.asarray(c) * np.asarray(w)
        )
        quad = Quadrature.gauss_legendre(8)
        space = _slab(quad)
        # anti-vacuity: |Ω·n| must genuinely vary on the face, else drop-abs
        # is a no-op.
        if np.ptp(np.abs(space.omega_dot_n[0])) <= 0.1:
            pytest.fail("dud config — |Ω·n| does not vary; abs-drop is a no-op")
        if np.array_equal(space.inner_product_weights, _inline_partial_current(space, quad)):
            pytest.fail("drop-abs did not bite — trace metric gate is vacuous")
        if not np.any(np.asarray(space.inner_product_weights) < 0.0):
            pytest.fail("drop-abs must yield a signed (negative-bearing) metric")


# ════════════════════════ Mode-11 execution sentinel ═══════════════════
class TestKernelActuallyExecuted:
    def test_trace_build_executes_the_kernel(self, monkeypatch):
        """The trace metric build routes through the module-global
        face_streaming_normal (a routed-around inline spelling would leave
        this sentinel unfired)."""
        calls: list[int] = []
        real = ats_mod.face_streaming_normal

        def wrapped(coeff, weight):
            calls.append(1)
            return real(coeff, weight)

        monkeypatch.setattr(ats_mod, "face_streaming_normal", wrapped)
        AngularTraceSpace.from_quadrature_and_layout(
            Quadrature.gauss_legendre(8),
            FaceLayout.from_named_shapes([("xmax", (8, 2))]),
        )
        if len(calls) == 0:
            pytest.fail(
                "trace metric build never called face_streaming_normal — the "
                "partial-current metric is not single-sourced through the kernel"
            )
