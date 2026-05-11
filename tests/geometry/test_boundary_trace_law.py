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

The TIGHTEST contract in this file is registry-isolation:
:class:`BoundaryTraceLaw` and the legacy
:class:`~orpheus.geometry.boundary.BoundaryOperator` MUST hold
SEPARATE class-level dicts. Without this guarantee, Wave 7's
rename of legacy concretes (``VacuumBoundaryOperator`` →
``VacuumInflow`` on the new hierarchy) would silently clobber the
legacy entries because :class:`~orpheus.numerics.registry.RegistryMixin`
keys by string, not by type.

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
def test_legacy_and_new_registries_are_disjoint() -> None:
    """The ``BoundaryOperator.registry`` (legacy) and
    ``BoundaryTraceLaw.registry`` (Wave 3) are separate dicts.

    This is the load-bearing assertion for Wave 7: when concrete
    BCs migrate onto :class:`BoundaryTraceLaw`, they MUST NOT
    leak into the legacy ``BoundaryOperator.registry`` (or vice
    versa), even though both hierarchies share
    :class:`~orpheus.numerics.registry.RegistryMixin` and could
    in principle clash on string keys (e.g. ``"vacuum"``).
    """
    # Distinct dict objects.
    assert (
        BoundaryOperator.registry is not BoundaryTraceLaw.registry
    )
    # Legacy registry still has its 6 concrete BCs.
    legacy_keys = set(BoundaryOperator.registry.keys())
    assert {
        "vacuum",
        "reflective",
        "white",
        "periodic",
        "albedo",
        "mixed",
    } <= legacy_keys
    # The new registry only has the test stub (and any future
    # subclasses we'll register in this test file). It does NOT
    # contain the legacy keys.
    new_keys = set(BoundaryTraceLaw.registry.keys())
    for legacy_key in ("vacuum", "reflective", "white", "periodic", "albedo", "mixed"):
        assert legacy_key not in new_keys
    # And the legacy registry has none of the new stub's key.
    assert "_stub_for_test" not in legacy_keys


@pytest.mark.foundation
def test_legacy_concrete_not_in_new_registry_values() -> None:
    """For maximum paranoia: no class registered against the
    legacy root is ALSO registered against the new root."""
    legacy_classes = set(BoundaryOperator.registry.values())
    new_classes = set(BoundaryTraceLaw.registry.values())
    assert legacy_classes.isdisjoint(new_classes)


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
