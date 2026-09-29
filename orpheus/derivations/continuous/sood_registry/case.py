r"""The schema of a published benchmark case: the configuration, its published values, and their citations.

:class:`La13511Case` holds one case of the registry: its cross sections, its
geometry kind and scattering order, the citation of the problem it poses, and
its published values (:class:`La13511Truth`) with their citations. The Sood,
Forster and Parsons cases (:mod:`.sood2003`) and the Atalay cases
(:mod:`.atalay1997`) are both built on it.

The two classes are transitional, and their names keep the 1999 report's
number until they retire. The reference-solution campaign
(``.claude/plans/reference_cache.md``, phase P4, the Sood step) splits a case
into the problem, a ``Specification``, and its answer, a ``PublishedSolution``;
the citations (:class:`~orpheus.data.citation.Citation`) already sit where
that design puts them: the problem's on the case, the values' on the truth.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Mapping

from orpheus.data.citation import Citation
from orpheus.data.macro_xs.mixture import Mixture
from orpheus.geometry import BC

if TYPE_CHECKING:
    from orpheus.geometry.structured_geometry import StructuredGeometry


@dataclass(frozen=True)
class La13511Truth:
    """Tabulated reference values for a Sood case.

    Different cases populate different subsets. ``None`` means the
    published reference does not tabulate that quantity.

    Parameters
    ----------
    k_eff_or_kinf : float
        Reference :math:`k_{\\rm eff}` (finite cases — usually 1.0,
        critical) or :math:`k_\\infty` (infinite cases).
    sources : tuple of Citation
        Where the values are printed, and the primary sources the
        publication credits for them: at least one. For a Sood case, the
        table or page of the 2003 paper, then the works its reference
        column names.
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
    critical_dimension_mfp : float | None
        Published critical dimension in mean free paths.

        For ``"slab"``: the half-thickness :math:`a` (F_N convention).
        For ``"sphere"`` / ``"cylinder"``: the radius :math:`R`.
        For ``"infinite"``: ``None``.

        Use case: registry-truth value (the published critical
        configuration). Multiply by :math:`1 / \\Sigma_t` to convert to
        cm; this is what :meth:`La13511Case.to_geometry` does internally.
        Living on Truth and not on geometry mirrors the architectural
        fact that this is a truth claim ("at this size, the configuration
        is critical"), not a geometric description.
    extrapolated_endpoint_mfp : float | None
        Published extrapolated endpoint :math:`z_0` in mean free paths.

        Standard transport-theory value used for diffusion-theory
        boundary conditions. Optional metadata; not all cases publish it.
    """

    k_eff_or_kinf: float
    sources: tuple[Citation, ...]
    flux_ratios: Mapping[float, float] | None = None
    flux_ratio_groupwise: Mapping[int, float] | None = None
    angular_flux_at_surface: Mapping[float, Mapping[float, float]] | None = None
    critical_dimension_mfp: float | None = None
    extrapolated_endpoint_mfp: float | None = None

    def __post_init__(self) -> None:
        if not self.sources:
            raise ValueError(
                "La13511Truth.sources is empty: a published value cites where "
                "it is printed"
            )


@dataclass(frozen=True)
class La13511Case:
    """A single Sood benchmark configuration.

    Production-protocol form: cross sections live in a
    :class:`Mixture` (the same object production solvers consume),
    geometry kind is a single string tag (``"slab"`` / ``"sphere"`` /
    ``"cylinder"`` / ``"infinite"``), and reference values live in a
    :class:`La13511Truth` (including the published critical dimension
    in mean free paths).

    Parameters
    ----------
    case_id : str
        The publication's identifier for the case (Sood's convention is
        ``<Material>-<Groups>-<Scattering>-<Geometry>``).
    problem : Citation
        The source that defines the problem, and where in it: for Sood,
        ``Citation("SoodForsterParsons2003", "problem N")``, N the problem
        number (1-75, the same in both editions).
    description : str
        One-line human-readable description.
    materials : dict[int, Mixture]
        Macroscopic cross sections keyed by material ID. Single-region
        cases use ``{0: Mixture(...)}``; multi-region cases (none in
        the first slice) add more keys.
    geometry_kind : str
        One of ``"slab"``, ``"sphere"``, ``"cylinder"``, ``"infinite"``.
        Use :meth:`to_geometry` to materialise a
        :class:`~orpheus.geometry.structured_geometry.StructuredGeometry`
        for finite kinds (raises for ``"infinite"``).
    scattering_order : int
        Legendre order of the scattering kernel (0 = isotropic,
        1 = P_1, 2 = P_2). All first-slice cases are isotropic
        (``= 0``).
    truth : La13511Truth
        All published reference values for this case (including
        :attr:`La13511Truth.critical_dimension_mfp` for finite cases).
    notes : str
        Remarks on the case: conventions, conversions, known issues.

    The published values' citations are on the truth
    (:attr:`La13511Truth.sources`): each citation goes where its claim
    lives. The free-text ``Provenance`` record this replaced retired in
    P1 step 2c.
    """

    case_id: str
    problem: Citation
    description: str
    materials: dict[int, Mixture]
    geometry_kind: str
    scattering_order: int
    truth: La13511Truth
    notes: str = ""

    def __post_init__(self) -> None:
        if self.geometry_kind not in {"slab", "sphere", "cylinder", "infinite"}:
            raise ValueError(
                f"La13511Case.geometry_kind must be one of "
                f"{{'slab', 'sphere', 'cylinder', 'infinite'}}; "
                f"got {self.geometry_kind!r}"
            )

    # ── New-API geometry adapter ───────────────────────────────────

    def to_geometry(self) -> "StructuredGeometry":
        r"""Convert this case's published mfp values to a cm-form
        :class:`StructuredGeometry` (the new geometry-layer object).

        Reads :attr:`La13511Truth.critical_dimension_mfp` and the
        primary material's :math:`\Sigma_t` (group 0) from
        ``self.materials[0].SigT[0]``. Computes the cm extent via
        :math:`\text{cm} = \text{mfp} / \Sigma_t`. Builds the coordinate
        system, the breakpoints of the one interval, and the boundary
        laws.

        Slab convention
        ---------------
        Returns the FULL slab width
        (:math:`2 \cdot \text{critical\_dimension\_mfp} / \Sigma_t`)
        with vacuum-vacuum BCs. This is the production-natural
        convention — F_N's natural half-thickness is recovered inside
        :class:`MomentSpace` from
        :attr:`StructuredGeometry.domain_extent_cm` ``/ 2``.

        Note: this encoding wastes the slab's natural half-symmetry —
        a future improvement would encode Sood symmetric slabs as
        half-slabs with reflective+vacuum BCs (half the cells in
        production solves at the same accuracy).

        Sphere / cylinder convention
        ----------------------------
        Single region of radius
        :math:`R = \text{critical\_dimension\_mfp} / \Sigma_t`,
        outer law vacuum; the centre of a solid cylinder or sphere is an
        interior point and carries no law, so the geometry takes the one
        outer law.

        Multi-region cases
        ------------------
        Every case of this registry is homogeneous, one material
        filling the body.

        Returns
        -------
        StructuredGeometry
            New-API geometry-layer object.

        Raises
        ------
        ValueError
            If :attr:`self.truth.critical_dimension_mfp` is None
            (e.g. for infinite-medium ``k_inf`` cases — those don't
            have a geometry; use :func:`solve_homogeneous_infinite` or
            :meth:`MomentSpace.solve_kinf` for those).
        ValueError
            If the case's geometry kind is not one of ``"slab"`` /
            ``"sphere"`` / ``"cylinder"`` (e.g. ``"infinite"``,
            ``"ISLC"``, or any future kind without a corresponding
            coordinate system).
        """
        # Local import to avoid a registry → geometry-layer import
        # cycle (the geometry layer doesn't know about cases).
        from orpheus.geometry import CoordSystem, StructuredGeometry

        kind = self.geometry_kind

        if kind == "infinite":
            raise ValueError(
                f"Case {self.case_id!r} is infinite-medium ("
                f"truth.k_eff_or_kinf = k_inf, no geometry). Use "
                f"orpheus.homogeneous.solver.solve_homogeneous_infinite "
                f"or MomentSpace.solve_kinf(mix) for k_inf cases — "
                f"to_geometry() is only defined for finite cases."
            )

        cd_mfp = self.truth.critical_dimension_mfp
        if cd_mfp is None:
            raise ValueError(
                f"Case {self.case_id!r}: truth.critical_dimension_mfp "
                f"is None — to_geometry() requires a published critical "
                f"dimension. (Multi-region cases without a single "
                f"published scalar are handled by sibling registry "
                f"builders, not by this method.)"
            )

        # cm = mfp / Σ_t (group 0 — the primary group's total XS sets
        # the mfp scale). For multi-group cases the choice of group is
        # by convention; Sood publishes mfp values in units of the
        # group-0 Σ_t throughout the catalogue.
        primary_mixture = next(iter(self.materials.values()))
        sigma_t = float(primary_mixture.SigT[0])
        cm = float(cd_mfp) / sigma_t

        # The registry's geometry kind → the geometry's coordinate system.
        _KIND_TO_COORD = {
            "slab": CoordSystem.CARTESIAN,
            "cylinder": CoordSystem.CYLINDRICAL,
            "sphere": CoordSystem.SPHERICAL,
        }
        if kind not in _KIND_TO_COORD:
            raise ValueError(
                f"Case {self.case_id!r}: geometry kind {kind!r} has no "
                f"coordinate system. Supported: {sorted(_KIND_TO_COORD)}."
            )
        coord = _KIND_TO_COORD[kind]

        # All first-slice cases are single-region with the primary
        # mixture at mat_id=0 (the convention this registry uses).
        mat_id = 0

        if coord is CoordSystem.CARTESIAN:
            # Full slab width (2 × half-thickness), vacuum-vacuum laws.
            return StructuredGeometry(
                coord=coord,
                breakpoints=(0.0, 2.0 * cm),
                mat_ids=(mat_id,),
                boundaries=(BC.vacuum, BC.vacuum),
            )

        # A solid cylinder or sphere of radius cm; the outer law vacuum.
        return StructuredGeometry(
            coord=coord,
            breakpoints=(0.0, cm),
            mat_ids=(mat_id,),
            boundaries=(BC.vacuum,),
        )


__all__ = ["La13511Case", "La13511Truth"]
