r"""The A|B|A cross-check problem and its harness, defined once (#405 P2 step 7b.2, spec §1.7b.2).

Fuel A | moderator B | fuel A at outer radii 0.5, 1.5, 2.0 cm, 2 groups,
reflective at r = R, on a sphere and on a cylinder: the problem of the
reference cross-check rows. It is posed as a
:class:`~orpheus.specification.specification.GeometrySpecification` from
ISOTROPIC (P0) mixtures: the reference solves isotropic scattering only (the
characteristic reference refuses an anisotropic emission; the trajectory
resolvent before it read P0 alone), and the SN rows run at
``scattering_order=0``, while the xs_library's ``get_mixture("B", "2g")``
carries a P1 moment (mean cosine 0.6). A specification built from those would
pose a problem neither side solves.

The harness of the rows that compare against the reference lives here too:
the reference at its working point (the characteristic reference since P1
step (d); its error estimates are
:mod:`tests.gates.derivations._characteristic_ladders`), the shape observables, the one
spelling of a tolerance scaled by a production reading, the cylinder's
verification and RECORD calls, and the strict-xfail mark.
"""
from __future__ import annotations

import functools
import math
from collections.abc import Mapping
from decimal import ROUND_DOWN, Decimal
from typing import TYPE_CHECKING, NamedTuple

import numpy as np
import pytest
import sympy

from orpheus.data.cells import CellCoefficient, Channel
from orpheus.data.macro_xs.mixture import Mixture
from orpheus.data.materials import Materials
from orpheus.derivations.common.xs_library import get_xs, make_mixture
from orpheus.geometry import BC, CoordSystem, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesh1D, Mesher
from orpheus.numerics.mesh_free_function import RegionwiseConstant, Symbolic
from orpheus.numerics.observable import Eigenvalue, FluxIntegral, Ratio
from orpheus.numerics.outcome import NotYet
from orpheus.numerics.question import Eigen
from orpheus.reference.solution import ReferenceSolution
from orpheus.reference.verification import ReferenceNotValid, verify_agreement
from orpheus.specification.specification import GeometrySpecification

if TYPE_CHECKING:
    from orpheus.derivations.continuous.characteristic import Resolution

#: The outer radius of each region, cm.
ABA_RADII = (0.5, 1.5, 2.0)
#: The material of each region: fuel A (id 0), moderator B (id 1), fuel A.
ABA_MATERIAL_IDS = (0, 1, 0)


def isotropic_mixture(key: str, groups: str = "2g") -> Mixture:
    """The xs_library's mixture ``key`` at ``groups`` with its P0 scattering only (no ``sig_s1``)."""
    xs = get_xs(key, groups)
    return make_mixture(sig_t=xs["sig_t"], sig_c=xs["sig_c"], sig_f=xs["sig_f"], nu=xs["nu"], chi=xs["chi"], sig_s=xs["sig_s"])


@functools.cache
def aba_materials() -> Materials:
    """Fuel A as material 0, moderator B as material 1, both isotropic."""
    return Materials({0: isotropic_mixture("A"), 1: isotropic_mixture("B")})


def aba_geometry(coord: CoordSystem) -> StructuredGeometry:
    """The A|B|A body on ``coord`` (spherical or cylindrical), reflective at r = R."""
    return StructuredGeometry.from_thicknesses(
        coord=coord, thicknesses=tuple(np.diff((0.0, *ABA_RADII))), mat_ids=ABA_MATERIAL_IDS,
        boundaries=(BC.reflective,),
    )


@functools.cache
def aba_specification(coord: CoordSystem) -> GeometrySpecification:
    """The k question on the A|B|A body: the specification the characteristic reference answers (P1 step (d))."""
    return GeometrySpecification(aba_materials(), aba_geometry(coord), Eigen(CellCoefficient.every(Channel.FISSION_EMISSION)))


def aba_uniform_width_mesh(coord: CoordSystem, n_cells: int) -> Mesh1D:
    """The A|B|A body meshed with ``n_cells`` cells of one width, each region holding its width's share of them."""
    geometry = aba_geometry(coord)
    outer = geometry.breakpoints[-1]
    return Mesher(geometry).partition(tuple(
        CellsByCount.uniform_width(round(n_cells * (b - a) / outer)) for a, b in geometry.intervals
    )).mesh


class AbaCrossSections(NamedTuple):
    """The A|B|A data per region (the direct-solver ladders' arrays): ``sigma_s`` is ``[region, from, to]``."""

    sigma_t: np.ndarray
    sigma_s: np.ndarray
    nu_sigma_f: np.ndarray
    chi: np.ndarray


def aba_xs_2g() -> AbaCrossSections:
    """Each region's isotropic data, read off :func:`aba_materials` (one definition of the data)."""
    mats = [aba_materials()[m] for m in ABA_MATERIAL_IDS]
    return AbaCrossSections(
        sigma_t=np.stack([np.asarray(m.SigT) for m in mats]),
        sigma_s=np.stack([m.SigS[0].toarray() for m in mats]),
        nu_sigma_f=np.stack([np.asarray(m.SigP) for m in mats]),
        chi=np.stack([np.asarray(m.chi) for m in mats]),
    )


# ── the reference ──────────────────────────────────────────────────────────


def aba_reference_at(coord: CoordSystem, resolution: Resolution) -> ReferenceSolution:
    """The characteristic reference on the A|B|A body at ``resolution``: lazy, uncertified (#566).

    Since P1 step (d) of ``.claude/plans/characteristic_reference_architecture.md`` (until then the trajectory
    resolvent, whose spelling the corroboration rows keep in ``tests/gates/derivations/_trajectory_resolvent_aba.py``).
    """
    from orpheus.derivations.continuous.characteristic import characteristic_reference

    return characteristic_reference(aba_specification(coord), resolution)


@functools.cache
def aba_reference(coord: CoordSystem) -> ReferenceSolution:
    """The reference at its working point (``_characteristic_ladders.WORKING_DEGREE``), once per session."""
    from tests.gates.derivations._characteristic_ladders import WORKING_DEGREE, rung

    return aba_reference_at(coord, rung(WORKING_DEGREE[_ABA_FIXTURE[coord]]))


#: Each body's fixture name in :mod:`tests.gates.derivations._characteristic_ladders`.
_ABA_FIXTURE = {CoordSystem.SPHERICAL: "aba_sphere", CoordSystem.CYLINDRICAL: "aba_cylinder"}


# ── the observables ────────────────────────────────────────────────────────


def fission_production_weight() -> RegionwiseConstant:
    """νΣ_f per region and group: the weight of the total fission production, the shape rows' gauge (step 8's
    ``Rate(fission production)`` resolved by hand, read through the mesh's region labels)."""
    return RegionwiseConstant(aba_xs_2g().nu_sigma_f)


def shape_observables(mesh: Mesh1D) -> tuple[tuple[int, int, Ratio], ...]:
    """Per SN cell i and group g: the cell average of φ_g over cell i gauged to unit total fission production,
    ``Ratio(FluxIntegral(indicator of cell i in group g / V_i), FluxIntegral(νΣ_f))``. Mesh-free observables
    (a weight of r), so the reference reads them in its own representation; together they are the shape
    metric of 2026-09-26 (fission-gauged cell averages over the SN cells), exactly ([M] probe5)."""
    n_groups = aba_specification(mesh.coord).n_groups
    denominator = FluxIntegral(fission_production_weight())
    edges, volumes = np.asarray(mesh.edges, dtype=float), np.asarray(mesh.volumes, dtype=float)
    r = Symbolic.r
    out = []
    for i in range(len(volumes)):
        a, b = sympy.Rational(float(edges[i])), sympy.Rational(float(edges[i + 1]))
        cell = sympy.Piecewise((sympy.Float(1.0 / float(volumes[i])), (r >= a) & (r < b)), (0, True))
        for g in range(n_groups):
            numerator = FluxIntegral(Symbolic.of(*(cell if h == g else 0 for h in range(n_groups))))
            out.append((i, g, Ratio(numerator, denominator)))
    return tuple(out)


# ── tolerances: a relative one, made absolute by a production reading ─────────


def truncated(value: float, figures: int = 3) -> float:
    """``value`` truncated toward zero to ``figures`` significant figures, never above it in magnitude: a scale that
    only TIGHTENS a tolerance it multiplies. The decimal truncation is taken on the value's shortest round-trip
    decimal (``repr``, so 1.4 stays 1.4) and the conversion back is stepped toward zero whenever it lands above the
    value (the floor-times-power-of-ten spelling it replaces rounded UP one ulp, ``truncated(1.4) =
    1.4000000000000001``, in 45 035 of 2e6 draws: qa F4)."""
    if value == 0.0 or not math.isfinite(value):
        raise ValueError(f"truncated: a scale is finite and non-zero, got {value!r}")
    exact = Decimal(repr(value))
    quantum = Decimal(1).scaleb(exact.copy_abs().adjusted() - figures + 1)
    result = float(exact.quantize(quantum, rounding=ROUND_DOWN))
    if abs(result) > abs(value):
        result = math.nextafter(result, 0.0)
    return result


def scaled_tolerance(relative: float, production_scale: float) -> float:
    """The absolute tolerance ``relative × truncated(production_scale)``: the verbs compare absolutely, and a
    cylinder xfail must refuse before the reference is read, so the scale is production's."""
    return relative * truncated(production_scale)


#: Production's algebraic error of a functional: no estimator exists yet (#564).
NO_ESTIMATOR = NotYet(564, "no algebraic-error estimator of a functional yet")


def verify_cylinder_k(solution, relative: float) -> None:
    """``verify_agreement`` of production's k against the cylinder reference, at ``relative`` scaled by production's
    k; the reference's refusal (no certificate) is the expected failure of every row that calls it."""
    k = solution.read(Eigenvalue()).value
    verify_agreement(
        solution, Eigenvalue(), aba_reference(CoordSystem.CYLINDRICAL), scaled_tolerance(relative, k), NO_ESTIMATOR,
    ).require()


# ── the RECORD ─────────────────────────────────────────────────────────────

#: RECORD values of the cylinder comparisons, k keys only (the user's cost ruling of 2026-10-03: the 80-ratio
#: shape set through the reference's extension costs about 38 min). ``[M]`` with
#: ``python -O -m tests.gates.derivations._characteristic_ladders records``: the SN keys 2026-09-26, unchanged since
#: (re-measured 2026-10-10: ``phase_c_k_sn`` within 1 ulp, ``unified_k`` 1.5e-10, ``standoff_sweep_k_nx40`` 0, all
#: inside the band); ``k_ref`` and ``phase_c_k_gap`` re-baselined 2026-10-10 at P1 step (d), when the reference
#: became the characteristic reference at the door default (until then the trajectory resolvent at (24, 16, 32):
#: ``k_ref`` 1.231036749830859, gap 5.97e-4; ``p1_step_d/ta_d/records.log``).
CYLINDER_3REG_RECORD: dict[str, float] = {
    "k_ref": 1.2317294528902538,
    "phase_c_k_sn": 1.2317720792844793,
    "phase_c_k_gap": 3.460574804554325e-05,
    "unified_k": 1.2310184196907399,
    "standoff_sweep_k_nx40": 1.23101841974857,
}
#: The eigenvalues are compared relatively, the gap absolutely.
CYLINDER_3REG_RECORD_RELATIVE = frozenset({"k_ref", "phase_c_k_sn", "unified_k", "standoff_sweep_k_nx40"})
#: Ten times the looser iterative floor: the 4x8 SN solves stop on a relative k increment of 1e-7 with a dominance
#: ratio near 0.95, so each k is within about 2e-6 of its fixed point.
CYLINDER_3REG_RECORD_BAND = 2e-5


def assert_record(readings: Mapping[str, float], recorded: Mapping[str, float], band: float, *, relative: frozenset[str]) -> None:
    """RECORD: each recorded reading still reads within ``band`` (relative for the names in ``relative``).

    Not verification: it pins what the code reads today, so a change on either side of a comparison that cannot
    yet verify reddens it. It raises (this module is a helper, so ``python -O`` would strip a bare ``assert``).
    """
    for name, value in recorded.items():
        scale = abs(value) if name in relative else 1.0
        if not abs(readings[name] - value) <= band * scale:
            raise AssertionError(
                f"record {name!r} moved: {readings[name]!r} against the recorded {value!r}. If the reference "
                f"changed (its working point, or the reference itself, #566), re-measure the records (python -O -m "
                f"tests.gates.derivations._characteristic_ladders records) and, once it is certified, "
                f"lift the xfails; otherwise an SN change moved the solve"
            )


def assert_cylinder_record(readings: Mapping[str, float]) -> None:
    """The cylinder RECORD over the keys ``readings`` carries, at its band and its relative set."""
    assert_record(
        readings, {name: CYLINDER_3REG_RECORD[name] for name in readings},
        CYLINDER_3REG_RECORD_BAND, relative=CYLINDER_3REG_RECORD_RELATIVE,
    )


#: The mark every row carries whose comparison needs the cylinder's certificate: strict (an XPASS fails) and
#: expecting only the verbs' own refusal, ``ReferenceNotValid`` (narrower than the ``AssertionError`` it had,
#: which any bare ``assert`` before the comparison also satisfied; ``vv-principles`` mode 8(4)).
awaits_cylinder_bound = pytest.mark.xfail(
    strict=True,
    raises=ReferenceNotValid,
    reason=(
        "#566 (#516 before P1 step (d)): the characteristic reference derives no bound, so it has no "
        "certificate and cannot anchor a verification; the verbs' refusal is the expected failure"
    ),
)
