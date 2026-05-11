r"""Foundation tests for the :class:`BoundaryTraceLaw` ABC.

This file pins the ABC contract shipped in Wave 3:

* The ABC is non-instantiable (``apply`` is abstract).
* A minimal concrete subclass with ``apply`` constructs and acts.
* The five ``assert_*`` invariants default to no-ops.
* The :pydata:`creates_sweep_cycle` ClassVar defaults to ``False``.
* :meth:`source` returns a :class:`NoSource` sentinel by default.
* :meth:`realize` raises :class:`NotImplementedError` (Wave 5 will
  wire the realiser registry).
* The :class:`BoundaryTraceLaw` registry self-populates via
  ``__init_subclass__(key=...)``.
* :attr:`geometry_map` and :attr:`response_kernel` default to
  ``None``.

Wave 7 update: the legacy ``BoundaryOperator`` ABC has been merged
into :class:`BoundaryTraceLaw`. The pre-Wave-7 registry-disjointness
tests were replaced with a UNIFIED-registry assertion — all 6 concrete
BCs (``vacuum`` / ``reflective`` / ``white`` / ``periodic`` /
``albedo`` / ``mixed``) now live in :pyattr:`BoundaryTraceLaw.registry`.
The legacy ``BoundaryOperator`` symbol is an alias of
:class:`BoundaryTraceLaw` (so ``BoundaryOperator.registry is
BoundaryTraceLaw.registry`` is now True).

Tagged ``@pytest.mark.foundation`` per :mod:`tests._harness`.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from orpheus.geometry.boundary import (
    BoundaryOperator,
    BoundarySource,
    BoundaryTraceLaw,
    ConstantInflowSource,
    NoSource,
)


# ─────────────────────────────────────────────────────────────────────
# Stub concrete law used across tests
# ─────────────────────────────────────────────────────────────────────


class _StubLaw(BoundaryTraceLaw, key="_stub_for_test"):
    """Minimal concrete law: scales the input by 0.5.

    The ``*args, **kwargs`` signature on ``apply`` matches the
    transition-period contract -- subclasses may carry the legacy
    2-arg form ``(psi_out, quadrature)`` during Waves 7-9 or the
    1-arg form ``(psi)`` post-Wave 10 with the same definition.
    """

    def apply(
        self, psi_out: np.ndarray, *args: Any, **kwargs: Any
    ) -> np.ndarray:
        return psi_out * 0.5


class _MockTrace:
    """Stand-in for :class:`InflowTraceSpace` with just a
    ``shape`` attribute -- enough to exercise
    :class:`NoSource` / :class:`ConstantInflowSource`."""

    def __init__(self, shape: tuple[int, ...]) -> None:
        self.shape = shape


# ─────────────────────────────────────────────────────────────────────
# ABC instantiability
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_abc_cannot_be_instantiated_directly() -> None:
    """The abstract ``apply`` blocks direct instantiation."""
    with pytest.raises(TypeError):
        BoundaryTraceLaw()  # type: ignore[abstract]


@pytest.mark.foundation
def test_minimal_concrete_subclass_works() -> None:
    """A subclass overriding ``apply`` constructs and acts."""
    law = _StubLaw()
    out = law.apply(np.ones(4))
    assert out.shape == (4,)
    np.testing.assert_array_equal(out, np.full(4, 0.5))


@pytest.mark.foundation
def test_concrete_subclass_callable_via_dunder() -> None:
    """``LinearOperatorMixin.__call__`` delegates to ``apply``."""
    law = _StubLaw()
    out = law(np.ones(4))
    np.testing.assert_array_equal(out, np.full(4, 0.5))


# ─────────────────────────────────────────────────────────────────────
# Default-property contract
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_geometry_map_default_is_none() -> None:
    assert _StubLaw().geometry_map is None


@pytest.mark.foundation
def test_response_kernel_default_is_none() -> None:
    assert _StubLaw().response_kernel is None


@pytest.mark.foundation
def test_source_default_is_no_source() -> None:
    law = _StubLaw()
    src = law.source
    assert isinstance(src, NoSource)
    trace = _MockTrace(shape=(8, 1))
    out = src.evaluate(trace)  # type: ignore[arg-type]
    assert out.shape == (8, 1)
    assert np.all(out == 0.0)


@pytest.mark.foundation
def test_creates_sweep_cycle_default_is_false() -> None:
    """Defaults to ``False``; reflective / periodic concretes in
    Wave 7 override to ``True``."""
    assert _StubLaw.creates_sweep_cycle is False
    assert BoundaryTraceLaw.creates_sweep_cycle is False


@pytest.mark.foundation
def test_capabilities_default_includes_apply() -> None:
    """The base ClassVar advertises ``CAP_APPLY``; concrete BCs
    may extend (e.g. specular adds ``CAP_APPLY_TRANSPOSE``)."""
    assert "apply" in _StubLaw.capabilities


@pytest.mark.foundation
def test_domain_and_range_default_to_none() -> None:
    """Defaults match :class:`LinearOperatorMixin` contract --
    operators predating function-space tagging return ``None``."""
    law = _StubLaw()
    assert law.domain is None
    assert law.range is None


# ─────────────────────────────────────────────────────────────────────
# assert_* invariants default to no-ops
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_assert_invariants_default_to_no_ops() -> None:
    """The five ``assert_*`` methods on the ABC default to no-ops.

    Concrete BCs in Wave 7 override the relevant ones with
    invariant logic that raises the typed errors from
    :mod:`orpheus.geometry.boundary._errors`. Default behaviour
    must not raise on ANY input.
    """
    law = _StubLaw()
    # Use ``None`` as the quadrature stand-in -- the default
    # implementations don't access it.
    assert law.assert_inflow_outflow_classification(None) is None  # type: ignore[arg-type]
    assert law.assert_outgoing_leakage_unconstrained(None) is None  # type: ignore[arg-type]
    assert law.assert_geometry_map_measure_preserving(None) is None  # type: ignore[arg-type]
    assert law.assert_response_positive_if_declared() is None
    assert law.assert_source_lives_on_incoming_trace(None) is None  # type: ignore[arg-type]


# ─────────────────────────────────────────────────────────────────────
# realize hook
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_realize_raises_not_implemented_in_wave_3() -> None:
    """The realiser registry ships in Wave 5; for Wave 3 the
    hook raises :class:`NotImplementedError` with a pointer to the
    plan."""
    law = _StubLaw()
    with pytest.raises(NotImplementedError, match="Wave 5"):
        law.realize(method_space=None)


# ─────────────────────────────────────────────────────────────────────
# Registry contracts
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_stub_law_self_registers_under_its_key() -> None:
    """The ``key="_stub_for_test"`` kwarg on the class statement
    populates :attr:`BoundaryTraceLaw.registry`."""
    assert BoundaryTraceLaw.registry.get("_stub_for_test") is _StubLaw
    # And we can construct via the registry's ``create`` classmethod.
    inst = BoundaryTraceLaw.create("_stub_for_test")
    assert isinstance(inst, _StubLaw)


@pytest.mark.foundation
def test_unified_registry_holds_all_concretes() -> None:
    """Post-Wave-7, the ONE registry on :class:`BoundaryTraceLaw`
    holds every concrete BC.

    Before Wave 7 the legacy ``BoundaryOperator`` ABC and
    :class:`BoundaryTraceLaw` held two disjoint registry dicts. Wave
    7 merged the ABCs (``BoundaryOperator`` is now an alias of
    :class:`BoundaryTraceLaw`), and the registry split with them —
    one dict containing the rank-1 concretes plus the test stub plus
    any Wave-7 additions (``prescribed_inflow``). Wave 11 removed
    ``MixedBoundaryOperator`` and its ``"mixed"`` registry key —
    rank-N compositions are now Wave-0 ``OperatorSum``-algebra over
    realised leaves, not a registered concrete BC.
    """
    # The two symbols are now the same class — registry IS the same dict.
    assert BoundaryOperator is BoundaryTraceLaw
    assert BoundaryOperator.registry is BoundaryTraceLaw.registry

    # Every Wave-4 rank-1 concrete BC lives in the unified registry.
    keys = set(BoundaryTraceLaw.registry.keys())
    expected_concretes = {
        "vacuum",
        "reflective",
        "white",
        "periodic",
        "albedo",
    }
    assert expected_concretes <= keys, (
        f"missing concretes in unified registry: "
        f"{expected_concretes - keys}"
    )
    # Wave 11 — ``"mixed"`` MUST no longer be a registered key.
    assert "mixed" not in keys, (
        "Wave 11 removed MixedBoundaryOperator; the 'mixed' key MUST "
        "no longer be in BoundaryTraceLaw.registry."
    )
    # The test stub registered in this module is also present.
    assert "_stub_for_test" in keys


# ─────────────────────────────────────────────────────────────────────
# Source helpers
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_no_source_evaluate_returns_zeros_of_trace_shape() -> None:
    src = NoSource()
    trace = _MockTrace(shape=(12, 2))
    out = src.evaluate(trace)  # type: ignore[arg-type]
    assert out.shape == (12, 2)
    assert out.dtype == np.float64
    assert np.all(out == 0.0)


@pytest.mark.foundation
def test_constant_inflow_source_returns_value_everywhere() -> None:
    src = ConstantInflowSource(value=2.5)
    trace = _MockTrace(shape=(6, 3))
    out = src.evaluate(trace)  # type: ignore[arg-type]
    assert out.shape == (6, 3)
    assert out.dtype == np.float64
    assert np.all(out == 2.5)


@pytest.mark.foundation
def test_boundary_source_protocol_is_runtime_checkable() -> None:
    """A bare instance of :class:`NoSource` or
    :class:`ConstantInflowSource` satisfies the Protocol via
    duck-typing (the ``evaluate`` method is present)."""
    assert isinstance(NoSource(), BoundarySource)
    assert isinstance(ConstantInflowSource(1.0), BoundarySource)
