"""The assembly functor's census over every consumer (numerics_untangling S3, gates G4b and G10).

Moved from ``tests/gates/numerics/test_assembly_functor.py`` by the layer DAG Phase 2 refresh (ruling 28). G10's iff
law runs here on the consumers' kinds and in the numerics file on the numerics kinds; "both kinds" counts their
union. The shared fixtures and the G10 law are :mod:`tests.gates.numerics._assembly_functor`."""

from __future__ import annotations

from typing import Any

import pytest

from orpheus.numerics.assembled_operator import SparseAssembledOperator, assemble
from orpheus.numerics.flat_operator import FlattenedOperator
from orpheus.numerics.operator import LinearOperator
from orpheus.transport.fields.cross_section_field import CrossSectionField
from tests.gates.numerics._assembly_functor import EXPECTED_ASSEMBLABLE, assert_is_assemblable_iff_assemble_succeeds, registry


pytestmark = [pytest.mark.foundation]


def _consumer_registry() -> dict[str, "Any"]:
    """The consumers' operator kinds (diffusion's parts and loss, the homogeneous loss and production), built
    lazily. The numerics kinds are :func:`tests.gates.numerics._assembly_functor.registry`; this census reads their union."""

    def diffusion():
        from tests.gates.diffusion.test_operators import _config_setup

        return _config_setup("albedo(0.3)/reflective")

    def diffusion_parts():
        from orpheus.diffusion.operators import DiffusionBoundaryOperator, LeakageOperator
        from orpheus.transport.operators.isotropic_transfer import IsotropicN2N, IsotropicScattering
        from orpheus.transport.operators.lift import BulkLift
        from orpheus.transport.operators.multiplication_operator import MultiplicationOperator

        mesh, mat_xs, template, _ = diffusion()
        ffs = mesh.full_field_space
        S = IsotropicScattering.from_material_xs(mat_xs, space=mesh.bulk_space)
        return {
            "L": LeakageOperator(mesh),
            "B": DiffusionBoundaryOperator(mesh),
            "C": MultiplicationOperator(CrossSectionField.total_of(mat_xs), domain=ffs, codomain=ffs),
            "S_plain": S,
            "N2N_plain": IsotropicN2N.from_material_xs(mat_xs, space=mesh.bulk_space),
            "S_lift": BulkLift(S, domain=ffs, codomain=ffs),
            "template": template,
            "mesh": mesh,
        }

    def homogeneous():
        from orpheus.derivations import get
        from orpheus.homogeneous.solver import HomogeneousProblem

        return HomogeneousProblem(next(iter(get("homo_2eg_n2n").materials.values())))

    return {
        "diffusion_L": lambda: diffusion_parts()["L"],
        "diffusion_B": lambda: diffusion_parts()["B"],
        "diffusion_C": lambda: diffusion_parts()["C"],
        "diffusion_S_plain": lambda: diffusion_parts()["S_plain"],
        "diffusion_N2N_plain": lambda: diffusion_parts()["N2N_plain"],
        "diffusion_S_lift": lambda: diffusion_parts()["S_lift"],
        "diffusion_loss_flat": lambda: (lambda p: FlattenedOperator(
            p["L"] + p["C"] - p["S_lift"] - p["B"], p["template"]))(diffusion_parts()),
        "homogeneous_loss": lambda: homogeneous().loss,
        "homogeneous_production": lambda: homogeneous().production,
    }


_CONSUMER_ASSEMBLABLE = {
    "diffusion_L", "diffusion_B", "diffusion_C", "diffusion_S_plain",
    "diffusion_N2N_plain", "diffusion_S_lift", "diffusion_loss_flat",
    "homogeneous_loss",
}


def _census_registry() -> dict[str, "Any"]:
    """The census's population: the numerics kinds and the consumers' kinds."""
    return {**registry(), **_consumer_registry()}


_CENSUS_ASSEMBLABLE = EXPECTED_ASSEMBLABLE | _CONSUMER_ASSEMBLABLE


class TestWhoConstructs:
    @pytest.mark.parametrize("subject", ["diffusion_loss", "homogeneous_loss"])
    def test_g4b_one_wrap_per_functor_call(self, subject: str, monkeypatch) -> None:
        """Transitively: a helper an emitter calls (``embed_bulk_assembly``,
        ``_bulk_assembly``) wrapping its result is invisible to G4a."""
        if subject == "diffusion_loss":
            from tests.gates.diffusion.test_operators import _config_setup, _loss

            mesh, mat_xs, template, _ = _config_setup("albedo(0.3)/reflective")
            op: LinearOperator = FlattenedOperator(_loss(mesh, mat_xs), template)
        else:
            from orpheus.derivations import get
            from orpheus.homogeneous.solver import HomogeneousProblem

            op = HomogeneousProblem(next(iter(get("homo_2eg_n2n").materials.values()))).loss
        wraps: list[int] = []
        original = SparseAssembledOperator.__init__

        def counting(self, *args, **kwargs):
            wraps.append(1)
            original(self, *args, **kwargs)

        monkeypatch.setattr(SparseAssembledOperator, "__init__", counting)
        assemble(op)
        assert len(wraps) == 1, f"{len(wraps)} carriers built for one assembly"


class TestTheFlagIsTheFunctorsDomain:
    def test_g10_the_registry_has_both_kinds(self) -> None:
        names = set(_census_registry())
        assert len(names) == 25
        assert _CENSUS_ASSEMBLABLE < names
        assert len(names - _CENSUS_ASSEMBLABLE) == 10

    @pytest.mark.parametrize("name", sorted(_consumer_registry()))
    def test_g10_is_assemblable_iff_assemble_succeeds(self, name: str) -> None:
        """G10 on this file's population; the law is stated once, in
        :func:`tests.gates.numerics._assembly_functor.assert_is_assemblable_iff_assemble_succeeds`."""
        assert_is_assemblable_iff_assemble_succeeds(name, _consumer_registry()[name](), _CENSUS_ASSEMBLABLE)


