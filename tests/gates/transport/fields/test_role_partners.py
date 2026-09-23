r"""The flux ↔ source/sink role pairs are a BIJECTION over the whole carrier
population (CS4c step 5, R-2 — G5.4 of the verification plan).

Every source/sink leaf declares its flux partner on its class statement
(``class AngularSourceSink(AngularField, flux=AngularFlux)``) and
:class:`~orpheus.transport.fields._bases.RolePair` registers BOTH
directions, so the operator tier can name the role it EMITS
(:meth:`into_role`) without an ``isinstance`` on the operand — the
carrier parse the step-5 census counted 12 times across three verbs.

**Two sources, two filters** (the 2026-08-19 countermeasure): the
declared partner map (the ``ClassVar`` the leaves set) against the
population derived by a runtime ``__subclasses__`` walk after importing
EVERY module of ``transport/{fields,source_sinks,residuals}`` explicitly.
⚠ `[M]` 2026-09-04: importing only the three package ``__init__``s
yields **19** leaves — ``AngularBoundaryFlux`` is not re-exported — so
the walk imports the modules, and that member is the positive control.

The population, stated (`[M]` 2026-09-04): **20 concrete leaves** — 7
flux, 7 source/sink, 5 residuals, 1 ``CrossSectionField``. The 7 pairs
are asserted in BOTH directions separately (a one-directional check is a
relation, not a bijection — vv anti-#14); the 6 unpaired leaves are
asserted to REFUSE, not silently dropped from the denominator.
"""
from __future__ import annotations

import importlib
import pkgutil
from dataclasses import dataclass
from typing import ClassVar

import numpy as np
import pytest

from orpheus.numerics.field import Field
from orpheus.numerics.space import FunctionSpace
from orpheus.numerics.units import SCALAR_FLUX_UNITS, SCALAR_RATE_UNITS, Unit
from orpheus.transport.fields._bases import FieldRole, RolePair, ScalarField

pytestmark = pytest.mark.foundation

_PACKAGES = (
    "orpheus.transport.fields",
    "orpheus.transport.source_sinks",
    "orpheus.transport.residuals",
)

#: The declared population — the second source the walk is checked against.
_FLUX = {
    "AngularFlux", "ScalarFlux", "HarmonicMomentFlux", "AngularBoundaryFlux",
    "ScalarBoundaryFlux", "RadialCharacteristicInteriorFlux",
    "RadialCharacteristicBoundaryFlux",
}
_SOURCE_SINK = {
    "AngularSourceSink", "ScalarSourceSink", "HarmonicMomentSourceSink",
    "AngularBoundarySourceSink", "ScalarBoundarySourceSink",
    "RadialCharacteristicInteriorSourceSink",
    "RadialCharacteristicBoundarySourceSink",
}
_UNPAIRED = {
    "AngularResidual", "ScalarResidual", "AngularBoundaryResidual",
    "RadialCharacteristicInteriorResidual",
    "RadialCharacteristicBoundaryResidual", "CrossSectionField",
}


def _all_subclasses(cls: type) -> set[type]:
    direct = set(cls.__subclasses__())
    return direct.union(*(_all_subclasses(s) for s in direct))


def _population() -> dict[str, type]:
    r"""Every CONCRETE transport leaf, by name, after importing every module."""
    for name in _PACKAGES:
        pkg = importlib.import_module(name)
        for module in pkgutil.iter_modules(pkg.__path__):
            importlib.import_module(f"{name}.{module.name}")
    leaves = {
        c.__name__: c
        for c in _all_subclasses(Field)
        if c.__module__.startswith("orpheus.transport.")
        and not c.__module__.endswith("._bases")
        and not c.__name__.startswith("_")
    }
    return leaves


def test_the_walk_sees_the_member_the_package_inits_drop():
    """Positive control: a package-``__init__``-only import reads 19 leaves."""
    assert "AngularBoundaryFlux" in _population()


def test_the_population_is_twenty_leaves_seven_of_each_role_six_unpaired():
    pop = _population()
    assert set(pop) == _FLUX | _SOURCE_SINK | _UNPAIRED, sorted(pop)
    assert len(pop) == 20


@pytest.mark.parametrize("sink_name", sorted(_SOURCE_SINK))
def test_every_source_sink_names_a_flux_and_the_flux_names_it_back(sink_name):
    pop = _population()
    sink = pop[sink_name]
    flux = sink.role_partner(FieldRole.FLUX)
    assert flux.__name__ in _FLUX, flux
    # the reverse direction, asserted separately (bijection, not relation)
    assert flux.role_partner(FieldRole.SOURCE_SINK) is sink
    # the identity halves
    assert sink.role_partner(FieldRole.SOURCE_SINK) is sink
    assert flux.role_partner(FieldRole.FLUX) is flux
    assert sink.role() is FieldRole.SOURCE_SINK
    assert flux.role() is FieldRole.FLUX


def test_the_seven_fluxes_are_hit_exactly_once():
    """Injectivity from the flux side: no two sinks share a flux."""
    pop = _population()
    hit = [pop[s].role_partner(FieldRole.FLUX).__name__ for s in sorted(_SOURCE_SINK)]
    assert sorted(hit) == sorted(_FLUX)


@pytest.mark.parametrize("name", sorted(_UNPAIRED))
def test_a_residual_or_coefficient_has_no_partner_and_says_so(name):
    cls = _population()[name]
    with pytest.raises(TypeError, match="no role partner"):
        cls.role_partner(FieldRole.FLUX)
    with pytest.raises(TypeError, match="no role"):
        cls.role()


def test_a_second_source_sink_for_one_flux_is_refused_at_class_creation():
    r"""The bijection is enforced at import time, not validated later."""
    from orpheus.transport.fields.scalar_flux import ScalarFlux

    with pytest.raises(TypeError, match="second one"):
        class _Twin(ScalarField, flux=ScalarFlux):  # noqa: F841
            UNITS: ClassVar[Unit] = SCALAR_RATE_UNITS


def test_flux_must_name_a_role_leaf():
    with pytest.raises(TypeError, match="flux= must name"):
        class _Bad(ScalarField, flux=int):  # type: ignore[arg-type]  # noqa: F841
            UNITS: ClassVar[Unit] = SCALAR_RATE_UNITS


class TestIntoRole:
    r"""``into_role`` — same space, same family fields, the other role's class."""

    def test_scalar_round_trip_carries_space_and_values(self):
        from orpheus.transport.fields.scalar_flux import ScalarFlux
        from orpheus.transport.source_sinks import ScalarSourceSink

        space = FunctionSpace(name="scalar", shape=(2, 3))
        phi = ScalarFlux(values=np.arange(6.0).reshape(2, 3), space=space)
        q = phi.into_role(FieldRole.SOURCE_SINK, phi.values * 2.0)
        assert type(q) is ScalarSourceSink
        assert q.space is phi.space
        np.testing.assert_array_equal(q.values, phi.values * 2.0)
        back = q.into_role(FieldRole.FLUX, q.values)
        assert type(back) is ScalarFlux and back.space is space
        same = q.into_role(FieldRole.SOURCE_SINK, q.values)
        assert type(same) is ScalarSourceSink

    def test_moment_family_fields_ride_across(self):
        r"""``L`` and ``spatial_moments`` are carried via ``dataclasses.fields``."""
        from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
        from orpheus.transport.source_sinks import HarmonicMomentSourceSink

        space = FunctionSpace(name="moments", shape=(2, 3, 2, 4))
        m = HarmonicMomentFlux(values=np.zeros((2, 3, 2, 4)), space=space, L=1)
        q = m.into_role(FieldRole.SOURCE_SINK, np.ones((2, 3, 2, 4)))
        assert type(q) is HarmonicMomentSourceSink
        assert q.L == 1 and q.spatial_moments == 1 and q.space is space

    def test_space_override_rebinds_the_output(self):
        from orpheus.transport.fields.scalar_flux import ScalarFlux

        a = FunctionSpace(name="scalar", shape=(2, 3))
        b = FunctionSpace(name="other", shape=(2, 3))
        phi = ScalarFlux(values=np.zeros((2, 3)), space=a)
        q = phi.into_role(FieldRole.SOURCE_SINK, np.zeros((2, 3)), space=b)
        assert q.space is b

    def test_an_unpaired_leaf_refuses(self):
        from orpheus.transport.fields.cross_section_field import CrossSectionField

        space = FunctionSpace(name="xs", shape=(2, 3))
        sig = CrossSectionField(values=np.ones((2, 3)), space=space)
        with pytest.raises(TypeError, match="no role partner"):
            sig.into_role(FieldRole.FLUX, sig.values)
