r"""Sood/Forster/Parsons LA-13511 (1999) benchmark case catalogue.

Production-protocol-aligned port of the legacy
``orpheus.derivations.continuous.fn_method.benchmarks.la13511`` module.
Each case now carries:

* **`materials: dict[int, Mixture]`** — keyed by integer material ID,
  exactly what :func:`orpheus.cp.solver.solve_cp` and
  :func:`orpheus.sn.solver.solve_sn` consume. Built via
  :func:`orpheus.derivations.common.xs_library.make_mixture` from the
  raw Sood XS components (ν, Σ_f, Σ_c, Σ_s, χ).
* **`mesh_template: MeshTemplate`** — encodes geometry kind +
  critical-dimension; ``mesh_template.build(n_cells)`` produces a
  :class:`orpheus.geometry.mesh.Mesh1D` at the desired refinement.
* **`truth: La13511Truth`** — all published reference values bundled
  in one struct. Different cases populate different subsets.

The legacy attributes (``case.sigma_t``, ``case.sigma_s`` etc.) are
retained as **read-only properties** that delegate to the materials
dict via :func:`.extractors.mixture_to_fn_arrays`. This keeps the F_N
test suite working unmodified while migration consumers flip to the
new shape one at a time.

ORPHEUS convention vs Sood convention
-------------------------------------

Sood et al. number energy groups :math:`g=N` (fast) → :math:`g=1`
(slow), the reverse of typical nuclear-engineering convention. This
module does the conversion at construction time; consumers see XS
arrays in **ORPHEUS convention** (``g=0`` fast, ``g=N-1`` slow).

The scattering matrix uses ORPHEUS's ``[from, to]`` convention:
``sigma_s[g, h]`` is :math:`\Sigma_{s, g \to h}`. The
``Sigma_g^{rem}`` removal cross sections used by Sood Eq 26-27 are
therefore ``sigma_t[g] - sigma_s[g, g]``.

Provenance
----------

* All XS values: LA-13511 Tables 1-67.
* All :math:`k_\infty`, :math:`r_c`, flux ratios: LA-13511 Section
  VII (1G), Section VII.B (2G), with primary references cited per case.

Python value precision: published values transcribed verbatim. Where
LA-13511 reports e.g. "k_inf = 2.612903" with 6 published digits, the
catalogue stores 2.612903 exactly as a Python float. Tolerance for
verification is 1e-5 (i.e., 5 significant figures match), so the
6th-digit rounding of the published value is not load-bearing.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping

import numpy as np

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.factories import Zone, mesh1d_from_zones
from orpheus.geometry.mesh import BC, Mesh1D


# ═══════════════════════════════════════════════════════════════════
# Mesh template + truth dataclasses
# ═══════════════════════════════════════════════════════════════════


_GEOMETRY_TO_COORD: dict[str, CoordSystem] = {
    "slab": CoordSystem.CARTESIAN,
    "sphere": CoordSystem.SPHERICAL,
    "cylinder": CoordSystem.CYLINDRICAL,
}


@dataclass(frozen=True)
class MeshTemplate:
    """Geometry + critical-dimension recipe for a Sood case.

    Method-agnostic: stores the shape of the domain (slab /
    sphere / cylinder / infinite / ISLC), the critical dimension as
    published, the group count, and the boundary conditions. The
    :meth:`build` method constructs a concrete :class:`Mesh1D` at the
    requested refinement using
    :func:`orpheus.geometry.factories.mesh1d_from_zones`.

    For ``geometry == "infinite"`` no mesh exists: :meth:`build`
    raises :class:`ValueError`. Infinite-medium consumers should use
    :func:`.builders.build_materials` and ignore the mesh entirely.

    Parameters
    ----------
    geometry : str
        One of ``"infinite"``, ``"slab"``, ``"sphere"``, ``"cylinder"``,
        ``"ISLC"``.
    critical_dimension_mfp : float | None
        Critical radius (sphere/cylinder) or half-thickness (slab) in
        mean free paths. ``None`` for infinite-medium cases.
    critical_dimension_cm : float | None
        Same as above but in cm.
    n_groups : int
        Number of energy groups (1, 2, 3, or 6).
    mat_id : int
        Material identifier used in the constructed mesh. The matching
        :class:`Mixture` lives in ``case.materials[mat_id]``.
    bc_left : BC
        Inner / left boundary condition. For sphere/cylinder default
        is :attr:`BC.reflective` (centreline); for slab default is
        :attr:`BC.vacuum`.
    bc_right : BC
        Outer / right boundary condition. Default
        :attr:`BC.vacuum`.
    """

    geometry: str
    critical_dimension_mfp: float | None
    critical_dimension_cm: float | None
    n_groups: int
    mat_id: int = 0
    bc_left: BC = field(default_factory=lambda: BC.vacuum)
    bc_right: BC = field(default_factory=lambda: BC.vacuum)

    def __post_init__(self) -> None:
        if self.geometry not in {"infinite", "slab", "sphere", "cylinder", "ISLC"}:
            raise ValueError(f"Unknown geometry {self.geometry!r}")
        if self.geometry == "infinite":
            if self.critical_dimension_cm is not None:
                raise ValueError(
                    "infinite geometry must have critical_dimension_cm=None"
                )
        else:
            if self.critical_dimension_cm is None:
                raise ValueError(
                    f"geometry={self.geometry!r} requires critical_dimension_cm"
                )

    @property
    def domain_extent_cm(self) -> float:
        """Total mesh extent in cm — what :meth:`build` actually constructs.

        Conventions:

        * **Slab**: ``2 * critical_dimension_cm`` (full symmetric slab
          ``[0, 2a]`` with vacuum BCs at both ends — F_N convention).
        * **Sphere / cylinder**: ``critical_dimension_cm`` (radius;
          mesh is ``[0, R]`` with reflective BC at the centre and
          vacuum at the outer surface unless overridden).
        * **Infinite / ISLC**: undefined, raises.
        """
        if self.geometry == "infinite":
            raise ValueError("infinite geometry has no domain_extent")
        if self.geometry == "ISLC":
            raise NotImplementedError("ISLC domain_extent not implemented")
        assert self.critical_dimension_cm is not None  # narrowed above
        if self.geometry == "slab":
            return 2.0 * float(self.critical_dimension_cm)
        return float(self.critical_dimension_cm)

    def build(self, n_cells: int = 64) -> Mesh1D:
        r"""Construct a :class:`Mesh1D` at the published critical dimension.

        Conventions (see :attr:`domain_extent_cm`):

        * **Slab** (``geometry == "slab"``): builds the FULL symmetric
          slab ``[0, 2a]`` where ``a = critical_dimension_cm``. This is
          the F_N method's natural domain — :math:`a` is the
          half-thickness, the published critical configuration is the
          full slab. Default BCs are vacuum at both ends.
        * **Sphere / cylinder**: builds ``[0, R]`` where
          ``R = critical_dimension_cm``. Default BCs: reflective at
          ``r = 0`` (centreline / axis), vacuum at the outer surface.

        Parameters
        ----------
        n_cells : int
            Number of equal-volume sub-cells.

        Returns
        -------
        Mesh1D
            A 1-D mesh with ``n_cells`` cells, homogeneous material
            ``mat_id``, and BCs from this template.

        Raises
        ------
        ValueError
            If ``geometry == "infinite"``.
        NotImplementedError
            If ``geometry == "ISLC"``.
        """
        if self.geometry == "infinite":
            raise ValueError(
                "infinite-medium cases have no mesh; consume materials directly"
            )
        if self.geometry == "ISLC":
            raise NotImplementedError(
                "ISLC mesh construction is not implemented in the first slice"
            )
        coord = _GEOMETRY_TO_COORD[self.geometry]
        zones = [Zone(
            outer_edge=self.domain_extent_cm,
            mat_id=self.mat_id,
            n_cells=n_cells,
        )]
        mesh = mesh1d_from_zones(zones, coord=coord)
        # Re-stamp with BC fields (mesh1d_from_zones doesn't set BCs).
        return Mesh1D(
            edges=mesh.edges,
            mat_ids=mesh.mat_ids,
            coord=mesh.coord,
            precomputed_volumes=mesh.precomputed_volumes,
            bc_left=self.bc_left,
            bc_right=self.bc_right,
        )


@dataclass(frozen=True)
class La13511Truth:
    """Tabulated reference values for a Sood LA-13511 case.

    Different cases populate different subsets. ``None`` means the
    published reference does not tabulate that quantity.

    Parameters
    ----------
    k_eff_or_kinf : float
        Reference :math:`k_{\\rm eff}` (finite cases — usually 1.0,
        critical) or :math:`k_\\infty` (infinite cases).
    flux_ratios : Mapping[float, float] | None
        For 1G cases with published flux table: dict mapping
        ``r/r_c`` to :math:`\\phi(r)/\\phi(0)`.
    flux_ratio_groupwise : Mapping[int, float] | None
        For 2G infinite cases: dict mapping ORPHEUS group index to
        :math:`\\phi_g/\\phi_0` (ratio relative to ORPHEUS group 0 =
        fast). For ``PU-2-0-IN``, this is
        ``{0: 1.0, 1: phi_slow/phi_fast}``.
    angular_flux_at_surface : Mapping[float, Mapping[float, float]] | None
        Reserved for future cases that publish surface angular flux
        :math:`\\psi(\\mu, r=R)`. ``None`` for first-slice cases.
    """

    k_eff_or_kinf: float
    flux_ratios: Mapping[float, float] | None = None
    flux_ratio_groupwise: Mapping[int, float] | None = None
    angular_flux_at_surface: Mapping[float, Mapping[float, float]] | None = None


# ═══════════════════════════════════════════════════════════════════
# Main case dataclass
# ═══════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class La13511Case:
    """A single LA-13511 benchmark configuration.

    Production-protocol form: cross sections live in a
    :class:`Mixture` (the same object production solvers consume),
    geometry lives in a :class:`MeshTemplate` (which builds a
    :class:`Mesh1D` on demand), reference values live in a
    :class:`La13511Truth`.

    Parameters
    ----------
    case_id : str
        Sood naming convention identifier
        (``<Material>-<Groups>-<Scattering>-<Geometry>``).
    problem_number : int
        Sequential problem number within LA-13511.
    description : str
        One-line human-readable description.
    materials : dict[int, Mixture]
        Macroscopic cross sections keyed by material ID. Single-region
        cases use ``{0: Mixture(...)}``; multi-region cases (none in
        the first slice) add more keys.
    mesh_template : MeshTemplate
        Geometry recipe + critical dimension. Use
        :meth:`MeshTemplate.build` to obtain a concrete mesh.
    scattering_order : int
        Legendre order of the scattering kernel (0 = isotropic,
        1 = P_1, 2 = P_2). All first-slice cases are isotropic
        (``= 0``).
    truth : La13511Truth
        All published reference values for this case.
    sood_table : int
        LA-13511 table number where this case is tabulated.
    primary_reference : str
        The peer-reviewed paper Sood cites as the source of the
        reference values.
    notes : str
        Any extra remarks (typo flags, conversion subtleties, etc).

    Notes
    -----
    The legacy F_N consumer interface (``case.sigma_t``,
    ``case.sigma_s``, ``case.nu_sigma_f``, ``case.chi``,
    ``case.geometry``, ``case.n_groups``,
    ``case.critical_dimension_mfp``, ``case.critical_dimension_cm``,
    ``case.k_eff_or_kinf``, ``case.flux_ratios``,
    ``case.flux_ratio_groupwise``) is exposed as **read-only
    properties** that delegate to ``materials[0]`` /
    ``mesh_template`` / ``truth``. This keeps the F_N test suite
    working unchanged.
    """

    case_id: str
    problem_number: int
    description: str
    materials: dict[int, Mixture]
    mesh_template: MeshTemplate
    scattering_order: int
    truth: La13511Truth
    sood_table: int
    primary_reference: str
    notes: str = ""

    # ── Legacy compatibility properties ────────────────────────────

    @property
    def n_groups(self) -> int:
        """Number of energy groups."""
        return self.mesh_template.n_groups

    @property
    def geometry(self) -> str:
        """Geometry kind (``"infinite"``, ``"slab"``, ...)."""
        return self.mesh_template.geometry

    @property
    def critical_dimension_mfp(self) -> float | None:
        """Critical dimension in mean free paths (or None for infinite)."""
        return self.mesh_template.critical_dimension_mfp

    @property
    def critical_dimension_cm(self) -> float | None:
        """Critical dimension in cm (or None for infinite)."""
        return self.mesh_template.critical_dimension_cm

    @property
    def k_eff_or_kinf(self) -> float:
        """Reference :math:`k_{\\rm eff}` (finite) or :math:`k_\\infty` (infinite)."""
        return self.truth.k_eff_or_kinf

    @property
    def flux_ratios(self) -> Mapping[float, float] | None:
        """Published flux ratio table, or None."""
        return self.truth.flux_ratios

    @property
    def flux_ratio_groupwise(self) -> Mapping[int, float] | None:
        """Published 2G groupwise flux ratio, or None."""
        return self.truth.flux_ratio_groupwise

    @property
    def sigma_t(self) -> np.ndarray:
        """Total XS in ORPHEUS convention, shape ``(n_groups,)``.

        Convenience accessor for legacy F_N consumers; pulls from
        ``self.materials[0]``.
        """
        from .extractors import mixture_to_fn_arrays
        sigma_t, _, _, _ = mixture_to_fn_arrays(self._primary_mixture)
        return sigma_t

    @property
    def sigma_s(self) -> np.ndarray:
        """P_0 scattering matrix in ``[from, to]`` convention, shape ``(n_groups, n_groups)``."""
        from .extractors import mixture_to_fn_arrays
        _, sigma_s, _, _ = mixture_to_fn_arrays(self._primary_mixture)
        return sigma_s

    @property
    def nu_sigma_f(self) -> np.ndarray:
        """Production XS :math:`\\nu \\Sigma_f`, shape ``(n_groups,)``."""
        from .extractors import mixture_to_fn_arrays
        _, _, nu_sigma_f, _ = mixture_to_fn_arrays(self._primary_mixture)
        return nu_sigma_f

    @property
    def chi(self) -> np.ndarray:
        """Fission spectrum, shape ``(n_groups,)``."""
        from .extractors import mixture_to_fn_arrays
        _, _, _, chi = mixture_to_fn_arrays(self._primary_mixture)
        return chi

    @property
    def _primary_mixture(self) -> Mixture:
        """The single mixture for first-slice (homogeneous) cases."""
        if len(self.materials) != 1:
            raise ValueError(
                f"Legacy property accessor expects single-region case; "
                f"case {self.case_id} has {len(self.materials)} materials. "
                f"Use case.materials[mat_id] explicitly."
            )
        return next(iter(self.materials.values()))


# ═══════════════════════════════════════════════════════════════════
# Helper: build a single-material 1G mixture from raw Sood XS
# ═══════════════════════════════════════════════════════════════════


def _mix_1g_isotropic(
    sigma_t: float,
    sigma_c: float,
    sigma_f: float,
    nu: float,
    sigma_s_self: float,
) -> Mixture:
    """Build a 1G isotropic Mixture from Sood's raw per-isotope XS components.

    Sood publishes :math:`(\\nu, \\Sigma_f, \\Sigma_c, \\Sigma_s, \\Sigma_t)`
    separately. Production solvers want the production XS
    :math:`\\nu \\Sigma_f` and the absorption XS components separately
    — :func:`make_mixture` accepts exactly that decomposition.
    """
    return make_mixture(
        sig_t=np.array([sigma_t]),
        sig_c=np.array([sigma_c]),
        sig_f=np.array([sigma_f]),
        nu=np.array([nu]),
        chi=np.array([1.0]),
        sig_s=np.array([[sigma_s_self]]),
    )


# ═══════════════════════════════════════════════════════════════════
# Case 1 — PUa-1-0-IN (Sood problem 1): 1G infinite medium, Pu-239 (a)
# ═══════════════════════════════════════════════════════════════════
#
# LA-13511 Table 2 (Pu-239 (a) cross sections, 1G isotropic):
#   ν = 3.24,  Σ_f = 0.0816,  Σ_c = 0.019584,  Σ_s = 0.225216,
#   Σ_t = 0.32640,  c = (Σ_s + νΣ_f)/Σ_t = 1.50.
# Reference: k_inf = 2.612903 (LA-13511 Eq 20).

PUA_1_0_IN = La13511Case(
    case_id="PUa-1-0-IN",
    problem_number=1,
    description="Pu-239 (a) bare infinite medium, 1G isotropic",
    materials={0: _mix_1g_isotropic(
        sigma_t=0.32640,
        sigma_c=0.019584,
        sigma_f=0.0816,
        nu=3.24,
        sigma_s_self=0.225216,
    )},
    mesh_template=MeshTemplate(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=1,
    ),
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=2.612903,
        flux_ratios=None,
    ),
    sood_table=2,
    primary_reference="LA-13511 Eq 20 (closed form)",
    notes=(
        "1G infinite-medium k_inf reduces to nu·Sigma_f/Sigma_a; the "
        "'c' factor in Eq 20 cancels algebraically (verified in "
        "fn_method.origins.k_inf_derivations.derive_kinf_1g_eq_20_simplifies_to_eq_19)."
    ),
)


# ═══════════════════════════════════════════════════════════════════
# Case 5 — PU-2-0-IN (Sood problem 44): 2G infinite medium, Pu
# ═══════════════════════════════════════════════════════════════════
#
# LA-13511 Tables 30-31 (Pu 2G isotropic, no upscatter):
#
#   Sood convention (g=2 fast, g=1 slow):
#     Fast (g=2):  ν_2=3.10,  Σ_2f=0.0936,  Σ_2c=0.00480,
#                  Σ_22s=0.0792, Σ_12s=0.0432 (from g=2 to g=1),
#                  Σ_2 = 0.2208,  χ_2 = 0.575.
#     Slow (g=1):  ν_1=2.93,  Σ_1f=0.08544, Σ_1c=0.0144,
#                  Σ_11s=0.23616, Σ_21s=0.0 (no upscatter),
#                  Σ_1 = 0.3360,  χ_1 = 0.425.
#
# In ORPHEUS convention (g=0 fast, g=1 slow):
#     sigma_t = [0.2208, 0.3360]
#     sigma_s[g=0=fast, g=0=fast] = Σ_22s = 0.0792
#     sigma_s[g=0=fast, g=1=slow] = (from fast to slow) = Σ_12s = 0.0432
#     sigma_s[g=1=slow, g=0=fast] = (from slow to fast) = Σ_21s = 0.0
#     sigma_s[g=1=slow, g=1=slow] = Σ_11s = 0.23616
#     sigma_f = [0.0936, 0.08544]
#     sigma_c = [0.00480, 0.0144]
#     nu = [3.10, 2.93]
#     chi = [0.575, 0.425]
#
# Reference: k_inf = 2.683767, φ_2/φ_1 = 0.675229 (Sood Eq 28-29 + Eq 32).

PU_2_0_IN = La13511Case(
    case_id="PU-2-0-IN",
    problem_number=44,
    description="Pu bare infinite medium, 2G isotropic, no upscatter",
    materials={0: make_mixture(
        sig_t=np.array([0.2208, 0.3360]),
        sig_c=np.array([0.00480, 0.0144]),
        sig_f=np.array([0.0936, 0.08544]),
        nu=np.array([3.10, 2.93]),
        chi=np.array([0.575, 0.425]),
        sig_s=np.array([
            [0.0792,  0.0432],   # from g=0 (fast):  → fast self, → slow downscatter
            [0.0,     0.23616],  # from g=1 (slow):  → fast upscatter (none), → slow self
        ]),
    )},
    mesh_template=MeshTemplate(
        geometry="infinite",
        critical_dimension_mfp=None,
        critical_dimension_cm=None,
        n_groups=2,
    ),
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=2.683767,
        flux_ratios=None,
        flux_ratio_groupwise={0: 1.0, 1: 1.0 / 0.675229},  # slow/fast = 1/Sood ratio
    ),
    sood_table=30,
    primary_reference="LA-13511 Eq 28-29 (k_inf) + Eq 32 (flux ratio)",
    notes=(
        "Sood Eq 28 has a typo: the chi_1 and chi_2 numerators have "
        "the wrong Sigma_g^rem factor — Eq 28 as printed reduces to "
        "2.862 (not 2.684) when Sigma_21s = 0, while Eq 29 (printed "
        "separately) reduces correctly. The SymPy derivation in "
        "fn_method.origins.k_inf_derivations.derive_kinf_2g_general_from_matrix "
        "exposes the correct general form by computing det(M)=0 from "
        "Eq 25 directly, and verifies the published Eq 29 against it."
    ),
)


# ═══════════════════════════════════════════════════════════════════
# Cases 2-4 — bare critical sphere/slab/cylinder of U-235 (a), 1G
# ═══════════════════════════════════════════════════════════════════
#
# Common XS (LA-13511 Tables 4-6, U-235 (a) 1G isotropic):
#   ν = 2.70,  Σ_f = 0.06528,  Σ_c = 0.013056,  Σ_s = 0.248064,
#   Σ_t = 0.32640,  c = (Σ_s + νΣ_f)/Σ_t = 1.30.
# All three cases are critical (k_eff = 1.0 by construction).

_UA_1G_KW = dict(
    sigma_t=0.32640,
    sigma_c=0.013056,
    sigma_f=0.06528,
    nu=2.70,
    sigma_s_self=0.248064,
)


# Case 2 — Ua-1-0-SL (Sood problem 12): 1G bare slab, U-235 (a)

UA_1_0_SL_STUB = La13511Case(
    case_id="Ua-1-0-SL",
    problem_number=12,
    description="U-235 (a) bare slab, 1G isotropic",
    materials={0: _mix_1g_isotropic(**_UA_1G_KW)},
    mesh_template=MeshTemplate(
        geometry="slab",
        critical_dimension_mfp=0.93772556,
        critical_dimension_cm=2.872934,
        n_groups=1,
        bc_left=BC.vacuum,
        bc_right=BC.vacuum,
    ),
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.0,
        flux_ratios={
            0.25: 0.9669506,
            0.50: 0.8686259,
            0.75: 0.7055218,
            1.00: 0.4461912,
        },
    ),
    sood_table=4,
    primary_reference="Kaper-Lindeman-Leaf 1974 NSE 54, 94",
    notes=(
        "Slab F_N solver shipped at ≤5e-6 absolute on a_c (see "
        "fn_method.slab.solve_fn_slab_bare_critical). Critical "
        "dimension is the half-thickness; slab full width is "
        "2*critical_dimension_cm = 5.745868 cm."
    ),
)


# Case 3 — Ua-1-0-CY (Sood problem 13): 1G bare cylinder, U-235 (a)

UA_1_0_CY_STUB = La13511Case(
    case_id="Ua-1-0-CY",
    problem_number=13,
    description="U-235 (a) bare infinite cylinder, 1G isotropic — STUB",
    materials={0: _mix_1g_isotropic(**_UA_1G_KW)},
    mesh_template=MeshTemplate(
        geometry="cylinder",
        critical_dimension_mfp=1.72500292,
        critical_dimension_cm=5.284935,
        n_groups=1,
        bc_left=BC.reflective,
        bc_right=BC.vacuum,
    ),
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.0,
    ),
    sood_table=5,
    primary_reference="Westfall-Metcalf 1973 NSE 52, 1",
    notes=(
        "STUB: F_N cylinder solver not yet implemented (needs "
        "Westfall-Metcalf 1973). Already cross-checked by Variant α "
        "at 8.5e-6."
    ),
)


# Case 4 — Ua-1-0-SP (Sood problem 14): 1G bare sphere, U-235 (a)

UA_1_0_SP_STUB = La13511Case(
    case_id="Ua-1-0-SP",
    problem_number=14,
    description="U-235 (a) bare sphere, 1G isotropic",
    materials={0: _mix_1g_isotropic(**_UA_1G_KW)},
    mesh_template=MeshTemplate(
        geometry="sphere",
        critical_dimension_mfp=2.4248249802,
        critical_dimension_cm=7.428998,
        n_groups=1,
        bc_left=BC.reflective,
        bc_right=BC.vacuum,
    ),
    scattering_order=0,
    truth=La13511Truth(
        k_eff_or_kinf=1.0,
    ),
    sood_table=6,
    primary_reference="Kaper-Lindeman-Leaf 1974 NSE 54, 94",
    notes=(
        "Sphere F_N solver shipped at ≤1e-7 absolute on R_c (see "
        "fn_method.sphere.solve_fn_sphere_bare_critical). Used as the "
        "structurally-independent L1 reference for Variant α sphere."
    ),
)


# ═══════════════════════════════════════════════════════════════════
# Public registry — name → case lookup
# ═══════════════════════════════════════════════════════════════════


ALL_FIRST_SLICE: tuple[La13511Case, ...] = (
    PUA_1_0_IN,
    PU_2_0_IN,
    UA_1_0_SL_STUB,
    UA_1_0_CY_STUB,
    UA_1_0_SP_STUB,
)
"""All first-slice cases. Cases 1+5 (k_inf) ship with full Branch-2
solvers; cases 2-4 are bare-critical slab/cylinder/sphere — slab and
sphere F_N shipped, cylinder F_N pending Westfall-Metcalf."""


LA13511_CASES: dict[str, La13511Case] = {
    case.case_id: case for case in ALL_FIRST_SLICE
}
"""Name → case mapping for ergonomic test access:

>>> from orpheus.derivations.continuous.sood_registry import LA13511_CASES
>>> case = LA13511_CASES["PUa-1-0-IN"]
>>> mixture = case.materials[0]
"""
