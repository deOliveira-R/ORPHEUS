r"""Solver adapters for the cross-method test protocol.

Each adapter wraps a continuous-reference solver to the
:class:`~tests.gates.cross_method.protocol.SolverAdapter` shape. The
adapter:

* reads cross sections directly off ``case.registry_case.materials``
  / ``case.materials`` (the production-protocol Mixture API);
* selects internal numerical parameters (n_modes for fn_method, the
  resolution rung for the characteristic reference) based on the
  requested case tolerance;
* performs unit conversions (mfp ↔ cm, half-thickness ↔ full
  slab);
* returns a :class:`ScalarResult` with the right ``tag``.

Phase D
-------

The pre-Phase-D ``mixture_to_fn_arrays`` extractor was retired as
part of the architectural reset; adapters now read
``mixture.SigT`` / ``SigS`` / ``SigP`` directly (the same pattern
the math-heart classes MomentSpace / Spectrum / BasisSpace
already use after their direct-__init__ migration).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


from .protocol import CrossMethodCase, ScalarResult, ScalarTag


# ═══════════════════════════════════════════════════════════════════
# F_N method adapters (fn_method package)
# ═══════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class FNSlabAdapter:
    r"""Adapter for :func:`...fn_method.slab.solve_fn_slab_bare_critical`.

    Reports the F_N method's predicted critical half-thickness in
    mean-free paths (the ``a_critical_mfp`` tag). Internally selects
    ``n_modes = 10`` by default — Grandjean-Siewert Table XI shows
    F_10 reaches ~5e-6 absolute on the slab half-thickness across
    the c-sweep, well below typical case tolerances.
    """

    name: str = "fn_slab"
    method: str = "fn_method"
    geometry: str = "slab"
    n_modes: int = 10

    def solve(self, case: CrossMethodCase) -> ScalarResult:
        from orpheus.derivations.continuous.fn_method.slab import (
            solve_fn_slab_bare_critical,
        )

        c = float(case.registry_case.materials[0].scattering_ratio[0])
        res = solve_fn_slab_bare_critical(c=c, n_modes=self.n_modes)
        return ScalarResult(
            tag="a_critical_mfp",
            value=float(res.a_critical_mfp),
            solver_name=self.name,
            metadata={
                "n_modes": self.n_modes,
                "determinant_residual": complex(res.determinant_residual),
                "nu0": float(res.nu0),
                "c": c,
            },
        )


@dataclass(frozen=True)
class FNSphereAdapter:
    r"""Adapter for :func:`...fn_method.sphere.solve_fn_sphere_bare_critical`.

    Reports ``R_critical_mfp``. Sphere F_N at ``n_modes = 10``
    reaches ~5e-8 absolute against Sood truth — exquisitely tight.
    """

    name: str = "fn_sphere"
    method: str = "fn_method"
    geometry: str = "sphere-1d"
    n_modes: int = 10

    def solve(self, case: CrossMethodCase) -> ScalarResult:
        from orpheus.derivations.continuous.fn_method.sphere import (
            solve_fn_sphere_bare_critical,
        )

        c = float(case.registry_case.materials[0].scattering_ratio[0])
        res = solve_fn_sphere_bare_critical(c=c, n_modes=self.n_modes)
        return ScalarResult(
            tag="R_critical_mfp",
            value=float(res.R_critical_mfp),
            solver_name=self.name,
            metadata={
                "n_modes": self.n_modes,
                "determinant_residual": complex(res.determinant_residual),
                "c": c,
            },
        )


@dataclass(frozen=True)
class FNReflectedSlabAdapter:
    r"""Adapter for :func:`...fn_method.slab.solve_fn_slab_reflected_critical`.

    Reflected-slab F_N (Neshat-Maiorino 1980). Returns
    ``tau_critical_mfp`` — the core half-thickness at criticality
    given the reflector configuration.

    Each case carries inline ``materials`` and a ``structured_geometry``
    that is a symmetric reflected slab (reflector, core, reflector). The
    adapter goes through :class:`MomentSpace`, the one door: it reads the
    body as a :class:`~orpheus.derivations.common.reference_body.ReflectedSlab`,
    requires one Sigma_t for both media, takes each ``c`` from
    :attr:`Mixture.scattering_ratio`, and converts the reflector width to
    mean free paths.

    There is currently no characteristic-reference counterpart for
    reflected slab — this adapter has no agreement partner. That
    one-sided coverage is intentional; see
    ``.claude/scratch/cross_method_test_protocol_assessment.md``
    §"Out of scope".
    """

    name: str = "fn_reflected_slab"
    method: str = "fn_method"
    geometry: str = "reflected-slab"
    n_modes: int = 7

    def solve(self, case: CrossMethodCase) -> ScalarResult:
        from orpheus.derivations.continuous.fn_method.moment_space import (
            MomentSpace,
        )

        # The reflected slab reaches the F_N solver through MomentSpace, the
        # one door (P1 step 2b): MomentSpace reads the geometry as a
        # symmetric reflected slab, checks that the core and the reflector
        # share one Sigma_t, and routes to solve_fn_slab_reflected_critical.
        if case.materials is None:
            raise ValueError(
                f"FNReflectedSlabAdapter: case {case.case_id!r} must carry "
                f"inline materials (core and reflector)."
            )
        solution = MomentSpace(
            geometry=_structured_geometry_for(case),
            materials=dict(case.materials),
            fn_order=self.n_modes,
        ).solve_critical()
        return ScalarResult(
            tag="tau_critical_mfp",
            value=float(solution.parameter_value),
            solver_name=self.name,
            metadata={
                "n_modes": self.n_modes,
                "c_core": solution.metadata["c_core"],
                "c_reflector": solution.metadata["c_reflector"],
                "reflector_half_thickness_mfp": solution.metadata["reflector_half_thickness_mfp"],
                "converged": bool(solution.converged),
            },
        )


# ═══════════════════════════════════════════════════════════════════
# Characteristic reference adapters (characteristic package)
# ═══════════════════════════════════════════════════════════════════


#: The working rung of the characteristic adapters: the joint ladder's panel degree
#: (:func:`tests.gates.derivations._characteristic_ladders.rung`). Every case's k at the rungs below and above it is
#: tabled in :mod:`.cases` (``CHARACTERISTIC_K``), where each case's tolerance is computed.
CHARACTERISTIC_WORKING_DEGREE = 4


@dataclass(frozen=True)
class CharacteristicAdapter:
    r"""Adapter for :func:`~orpheus.derivations.continuous.characteristic.characteristic_reference`: k of the case's
    own geometry and materials, posed as a specification.

    The successor of the three trajectory-resolvent adapters (P1 step (e1b) of
    ``.claude/plans/characteristic_reference_architecture.md``, the user's ruling 3 of 2026-10-10). One class serves
    every geometry: the reference reads the chart and the walls off the case's :class:`StructuredGeometry`, so the
    slab, the sphere and the closed sphere differ only in the name the cases key their tolerances on and in the
    tag (``k_inf`` for the closed body, whose k is the medium's). No convention is converted here: the geometry is
    in centimetres and the published critical dimension is read off it by the cases.
    """

    name: str
    geometry: str
    tag: ScalarTag
    method: str = "characteristic"
    degree: int = CHARACTERISTIC_WORKING_DEGREE

    def solve(self, case: CrossMethodCase) -> ScalarResult:
        from orpheus.data.cells import CellCoefficient, Channel
        from orpheus.data.materials import Materials
        from orpheus.derivations.continuous.characteristic import characteristic_reference
        from orpheus.numerics.observable import Eigenvalue
        from orpheus.numerics.question import Eigen
        from orpheus.specification.specification import GeometrySpecification
        from tests.gates.derivations._characteristic_ladders import rung

        if case.materials is not None:
            materials = case.materials
        elif case.registry_case is not None:
            materials = case.registry_case.materials
        else:
            raise ValueError(f"CrossMethodCase {case.case_id!r} carries neither inline materials nor a registry case")
        specification = GeometrySpecification(
            Materials(dict(materials)), _structured_geometry_for(case), Eigen(CellCoefficient.every(Channel.FISSION_EMISSION)),
        )
        k = characteristic_reference(specification, rung(self.degree)).read(Eigenvalue()).value
        return ScalarResult(tag=self.tag, value=float(k), solver_name=self.name, metadata={"degree": self.degree})


CHARACTERISTIC_SLAB = CharacteristicAdapter("characteristic_slab", "slab", "k_eff")
CHARACTERISTIC_SPHERE = CharacteristicAdapter("characteristic_sphere", "sphere-1d", "k_eff")
CHARACTERISTIC_SPHERE_CLOSED = CharacteristicAdapter("characteristic_sphere_closed", "closed-sphere-1d", "k_inf")


# ═══════════════════════════════════════════════════════════════════
# Helpers — XS / parameter extraction from CrossMethodCase
# ═══════════════════════════════════════════════════════════════════


def _extract_1g_xs(case: CrossMethodCase) -> tuple[float, float, float]:
    r"""Extract :math:`(\Sigma_t, \Sigma_s, \nu\Sigma_f)` for a 1G case
    from a registry-backed case.

    Pulls from ``case.registry_case.materials[0]`` via
    :func:`mixture_to_fn_arrays`. Raises if the case is multi-group
    (1G adapters can't consume those) or if the case carries no
    registry case.
    """
    if case.registry_case is None:
        raise ValueError(
            f"CrossMethodCase {case.case_id!r} has registry_case=None; "
            f"the registry-backed XS extractor cannot serve this case."
        )
    return _xs_from_materials_dict(
        case.registry_case.materials, case.case_id
    )


def _xs_from_materials_dict(
    materials: dict, case_id: str,
) -> tuple[float, float, float]:
    """Common backend: pull 1G ``(σ_t, σ_s, νσ_f)`` from a materials dict.

    Reads directly off ``Mixture.SigT`` / ``SigS[0]`` / ``SigP``
    (the production-protocol surface).
    """
    primary = materials[0]
    sigma_t_arr = np.asarray(primary.SigT, dtype=float)
    sigma_s_arr = primary.SigS[0].toarray().astype(float)
    nu_sigma_f_arr = np.asarray(primary.SigP, dtype=float)
    if sigma_t_arr.shape[0] != 1:
        raise ValueError(
            f"_xs_from_materials_dict: case {case_id!r} is "
            f"{sigma_t_arr.shape[0]}G; expected 1G"
        )
    return (
        float(sigma_t_arr[0]),
        float(sigma_s_arr[0, 0]),
        float(nu_sigma_f_arr[0]),
    )


def _structured_geometry_for(case: CrossMethodCase):
    r"""Resolve the :class:`StructuredGeometry` for a case.

    Reads from ``case.structured_geometry`` (inline / override path)
    or builds it via ``case.registry_case.to_geometry()`` (registry
    path), whichever is populated. The override path takes precedence
    when both are present (cross-method agreement gates substitute a
    predicted critical dimension by setting an inline structured
    geometry).
    """
    if case.structured_geometry is not None:
        return case.structured_geometry
    if case.registry_case is not None and hasattr(
        case.registry_case, "to_geometry"
    ):
        return case.registry_case.to_geometry()
    raise ValueError(
        f"_structured_geometry_for: case {case.case_id!r} has neither "
        f"inline structured_geometry nor a registry_case carrying one."
    )


# ═══════════════════════════════════════════════════════════════════
# Adapter registry — used by tests and (future) agreement-matrix renderer
# ═══════════════════════════════════════════════════════════════════


ADAPTERS_BY_NAME: dict[str, object] = {
    "fn_slab": FNSlabAdapter(),
    "fn_sphere": FNSphereAdapter(),
    "fn_reflected_slab": FNReflectedSlabAdapter(),
    "characteristic_slab": CHARACTERISTIC_SLAB,
    "characteristic_sphere": CHARACTERISTIC_SPHERE,
    "characteristic_sphere_closed": CHARACTERISTIC_SPHERE_CLOSED,
}
"""All registered adapters. New adapters MUST register here so the
agreement-matrix renderer (and future cross-method audit tools)
can discover them.
"""
