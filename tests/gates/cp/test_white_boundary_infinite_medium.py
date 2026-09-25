r"""A homogeneous CP cell with a white boundary is an infinite medium.

The closed form (Wigner–Seitz). A homogeneous cell whose boundary
re-emits its outgoing partial current isotropically (the white, or Mark,
boundary) cannot lose neutrons, so it is an infinite medium: in every
group the flux is flat in space, the group spectrum is the
infinite-medium eigenvector, and ``k = k∞``. For discrete collision
probabilities the statement holds EXACTLY, not to discretisation error:
the white fill makes every row of :math:`P^\infty` sum to 1
(:eq:`p-inf`, :eq:`cp-infinite-lattice-sum`), reciprocity
:math:`\Sigma_t V_i P_{ij} = \Sigma_t V_j P_{ji}` then makes the uniform
vector a fixed point of the CP equations, and the solve reduces to the
0-D eigenproblem. So the rows below compare a CP solve with a closed
form to the solver's own tolerance.

The reference. The dense 0-D pencil ``ρ(A⁻¹F)`` with
``A = diag(Σ_t) − Σ_s0ᵀ`` (``SigS`` is stored ``[from, to]``) and
``F = χ ⊗ νΣ_f``, assembled here from the raw arrays of the
:class:`~orpheus.data.macro_xs.mixture.Mixture`, independently of every CP
matrix. The dense eigensolver is a shared, trusted primitive
(``instrument-doctrine`` X4); the pencil itself is pinned against the
SymPy k∞ of the homogeneous derivations by
:func:`test_the_pencil_reproduces_the_sympy_kinf`.

The rows, per configuration:

- C1, ``|k − k∞| ≤ 10 × keff_tol``, with ``keff_tol`` read from the
  :class:`~orpheus.cp.solver.CPParams` that drove the solve.
- C2, per group, ``max_i |φ_ig / ⟨φ_g⟩ − 1| ≤ _FLATNESS_BAND``, with
  ``⟨φ_g⟩`` the volume-weighted cell average.
- C3, the group spectrum ``⟨φ_g⟩ / Σ_g ⟨φ_g⟩`` equals the pencil's
  eigenvector to ``10 × flux_tol``.

The three rows partition the defects, measured 2026-09-25 with four
in-process arms over this module (``python -O``; per-arm counts are over
the 33 non-#509 configurations, 27 for C2):

========================================  ======  ======  ======
arm                                       C1      C2      C3
========================================  ======  ======  ======
A: the white fill dropped (vacuum)        33 red  27 red  33 red
B: the white fill halved (albedo ½)       33 red  27 red  33 red
C: a 1e-3 row-sum defect, one P_cell row  blind   27 red  blind
D: the scattering transpose flipped       33 red  blind   33 red
========================================  ======  ======  ======

So ``k`` and the spectrum are invariant under a conservation defect of
one row of :math:`P^{cell}` (``instrument-doctrine``, the stabiliser;
``vv-principles`` Mode 12): C2 is the only row that sees it, and C1
alone would be a gate whose invariance group contains that error. The
flat shape is invariant under a wrong scattering convention, which only
C1 and C3 see.

What they null: spatial heterogeneity (one material). The heterogeneous
continuous rung is #507; the flat-source CP discrete-exact cases carry
the heterogeneous discrete claim.

The grid: coordinate ∈ {slab, cylinder, sphere} × R ∈ {0.5, 2, 10} cm
(mixture A has ``Σ_t`` between 0.4 and 1.05 cm⁻¹, so these are about
0.2 to 10 mean free paths; 0.5 is leakage-dominated, where the white
boundary carries the whole answer) × groups ∈ {2, 4}. At one group
``k = νΣ_f/Σ_a`` is independent of the flux shape (``vv-principles``
anti-pattern 3), and the 1G sphere at R = 10 shows ``k = k∞`` exactly
while its shape is wrong (the RECORD row below). Also: a geometrically
graded mesh (11 cells, refined toward the OUTER boundary) at R = 0.5 and
R = 4, and a ONE-cell mesh at R = 2. A one-cell mesh is flat by
construction, so it carries C1 and C3 and not C2.

#509. The CP sphere with a white boundary is not flat at large optical
radius, and at two or more groups its k and spectrum are wrong. Measured
2026-09-25 (``python -O``, ``keff_tol = flux_tol = 1e-12``): green at the
largest group optical radius ``max_g Σ_t R`` = 4.0 (2G, graded, R = 4);
red at 4.2 (4G, graded, R = 4: ``k/k∞ − 1 = 4.4e-3``, flatness 2.9e-3)
and at 10 and 10.5 (R = 10: ``k/k∞ − 1 = 4.7e-2`` at 2G, ``6.9e-2`` at
4G). Those sphere cells are strict xfails citing #509;
:func:`test_the_509_cells_are_strict_xfails` asserts the strictness by
introspection, and :func:`test_sphere_R10_flatness_defect_is_present`
records the defect's magnitude so that only a full fix turns it red.
"""
from __future__ import annotations

import contextlib
import dataclasses
import functools
import io

import numpy as np
import pytest

from orpheus.cp.solver import CPParams, CPResult, solve_cp
from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations import get as get_case
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, CoordSystem, Mesh1D

_MAT_ID = 2

#: The solver tolerances that drive every solve here. C1 and C3 read
#: their bounds off the instance that drove the solve, never a copy.
_PARAMS = CPParams(keff_tol=1e-12, flux_tol=1e-12, max_outer=2000)

#: C2's band: a relative flatness bound measured over the population, not
#: derived from one tolerance (the flux is flat to rounding, independently
#: of how far the power iteration ran). Measured 2026-09-25 over the 27
#: spatial configurations below that realise the closed form: worst
#: 5.4e-14 (CAR-R0.5-graded11-4g), so the band has about 19x headroom.
#: The weakest signal it must see is 3.4e-5 (arm C, a 1e-3 row-sum defect
#: on one row of P_cell, weakest over the same 27 configurations).
_FLATNESS_BAND = 1e-12

_509_REASON = (
    "#509: CP sphere white boundary not flat at group optical radius "
    "max_g Sigma_t R >= about 4.2 (measured 2026-09-25)"
)

_FOUNDATION = "tests/gates/cp/test_verification.py::TestMultiGroupProperties"
_PENCIL_ROWS = (
    "tests/gates/cp/test_white_boundary_infinite_medium.py"
    "::test_the_pencil_reproduces_the_sympy_kinf",
)


@dataclasses.dataclass(frozen=True)
class WhiteCell:
    """One homogeneous white-boundary configuration of mixture A."""

    coord: CoordSystem
    radius: float  # outer radius or slab width [cm]
    mesh: str  # "uniform", "graded" or "one-cell"
    n_cells: int
    groups: str  # xs_library group key, "2g" or "4g"

    @property
    def id(self) -> str:
        return (f"{self.coord.name[:3]}-R{self.radius:g}-{self.mesh}"
                f"{self.n_cells}-{self.groups}")

    def edges(self) -> np.ndarray:
        """Cell edges [cm]; the graded mesh is geometric, refined toward
        the outer boundary (spacing ratio about 0.66 per cell)."""
        if self.mesh == "uniform":
            return np.linspace(0.0, self.radius, self.n_cells + 1)
        if self.mesh == "graded":
            edges = self.radius * (1 - np.geomspace(1.0, 1e-2, self.n_cells + 1)) / (1 - 1e-2)
            edges[0], edges[-1] = 0.0, self.radius
            return edges
        return np.array([0.0, self.radius])

    def mesh1d(self) -> Mesh1D:
        edges = self.edges()
        # The slab has two faces; the curvilinear cell one (its axis or
        # centre is not a boundary).
        bc_left = BC.white if self.coord is CoordSystem.CARTESIAN else None
        return Mesh1D(edges=edges, mat_ids=np.full(len(edges) - 1, _MAT_ID),
                      coord=self.coord, bc_left=bc_left, bc_right=BC.white)

    def foundation_rungs(self) -> tuple[str, ...]:
        """The CP matrix-property rungs for this coordinate (2G): the
        white-filled P_inf row sums and a reciprocity pair."""
        tag = f"2g-CoordSystem.{self.coord.name}"
        return (f"{_FOUNDATION}::test_row_sums_multigroup[{tag}]",
                f"{_FOUNDATION}::test_reciprocity_multigroup[{tag}]")


def _grid() -> list[WhiteCell]:
    cells = []
    for coord in (CoordSystem.CARTESIAN, CoordSystem.CYLINDRICAL, CoordSystem.SPHERICAL):
        for groups in ("2g", "4g"):
            cells += [WhiteCell(coord, radius, "uniform", n, groups)
                      for radius, n in ((0.5, 5), (2.0, 10), (10.0, 20))]
            cells += [WhiteCell(coord, radius, "graded", 11, groups) for radius in (0.5, 4.0)]
            cells.append(WhiteCell(coord, 2.0, "one-cell", 1, groups))
    return cells


#: The sphere cells measured red on today's tree (#509), by id. The group
#: optical radius max_g Σ_t R is 10, 10.5 and 4.2 for these; the largest
#: sphere cell measured green is 4.0 (SPH-R4-graded11-2g).
_509_CELLS = frozenset({
    "SPH-R10-uniform20-2g",
    "SPH-R10-uniform20-4g",
    "SPH-R4-graded11-4g",
})


def _is_509(cell: WhiteCell) -> bool:
    return cell.id in _509_CELLS


def _params(cells: list[WhiteCell]) -> list:
    out = []
    for cell in cells:
        marks = [pytest.mark.rests_on(*cell.foundation_rungs(), *_PENCIL_ROWS)]
        if _is_509(cell):
            marks.append(pytest.mark.xfail(strict=True, reason=_509_REASON))
        out.append(pytest.param(cell, id=cell.id, marks=marks))
    return out


@functools.cache
def _mixture(groups: str) -> Mixture:
    return get_mixture("A", groups)


_ALL_CELLS = _grid()
_SPATIAL_CELLS = [c for c in _ALL_CELLS if c.mesh != "one-cell"]
_ALL_PARAMS = _params(_ALL_CELLS)
_SPATIAL_PARAMS = _params(_SPATIAL_CELLS)


@functools.cache
def _infinite_medium(groups: str) -> tuple[float, np.ndarray]:
    """k∞ and the normalised group spectrum of mixture A, from the dense
    0-D pencil assembled from the raw arrays."""
    mix = _mixture(groups)
    for order, sig2 in enumerate(mix.Sig2):
        # Premise: A below omits the (n,2n) loss; every library mixture
        # ships Sig2 = 0, and a non-zero one would make the pencil wrong.
        assert sig2.count_nonzero() == 0, (
            f"mixture A {groups} carries (n,2n) at order {order}; the pencil omits it"
        )
    sig_s0 = mix.SigS[0].toarray()  # [from, to] [1/cm]
    loss = np.diag(mix.SigT) - sig_s0.T  # [1/cm]
    production = np.outer(np.asarray(mix.chi), mix.SigP)  # [1/cm]
    eigenvalues, eigenvectors = np.linalg.eig(np.linalg.solve(loss, production))
    dominant = int(np.argmax(np.real(eigenvalues)))
    spectrum = np.abs(np.real(eigenvectors[:, dominant]))
    return float(np.real(eigenvalues[dominant])), spectrum / spectrum.sum()


@functools.cache
def _solve(cell: WhiteCell) -> tuple[CPResult, CPParams]:
    params = _PARAMS
    with contextlib.redirect_stdout(io.StringIO()):
        result = solve_cp({_MAT_ID: _mixture(cell.groups)}, cell.mesh1d(), params)
    assert result.record.fully_converged, f"{cell.id}: the CP solve did not converge"
    return result, params


def _group_averages(cell: WhiteCell, result: CPResult) -> np.ndarray:
    """The volume-weighted cell average of the flux per group, ``(ng,)``."""
    volumes = result.geometry.volumes
    return volumes @ result.flux / volumes.sum()


@pytest.mark.foundation
@pytest.mark.rests_on(
    "tests/gates/homogeneous/test_homogeneous.py::test_kinf_exact[homo_2eg]",
    "tests/gates/homogeneous/test_homogeneous.py::test_kinf_exact[homo_4eg]",
)
@pytest.mark.parametrize("groups, case_name", [("2g", "homo_2eg"), ("4g", "homo_4eg")])
def test_the_pencil_reproduces_the_sympy_kinf(groups: str, case_name: str) -> None:
    r"""The in-test reference agrees with the SymPy k∞ of the homogeneous
    derivations for the same mixture (a different derivation: the
    characteristic polynomial, not a dense eigensolve), so the rows below
    compare against a reference that is itself anchored. The case's
    material is checked to be mixture A first."""
    case = get_case(case_name)
    (case_mix,) = case.materials.values()
    mix = _mixture(groups)
    for field in ("SigT", "SigP"):
        np.testing.assert_array_equal(getattr(case_mix, field), getattr(mix, field))
    k_inf, _ = _infinite_medium(groups)
    assert abs(k_inf - case.k_inf) <= 1e-13 * case.k_inf, (k_inf, case.k_inf)


@pytest.mark.l1
@pytest.mark.verifies("cp-white-cell-infinite-medium", "p-inf", "cp-infinite-lattice-sum")
@pytest.mark.parametrize("cell", _ALL_PARAMS)
def test_k_is_the_infinite_medium_k(cell: WhiteCell) -> None:
    """C1: ``|k − k∞| ≤ 10 × keff_tol`` (keff_tol of the driving
    CPParams). Measured worst 2.1e-13 (CAR-R4-graded11-2g) against 1e-11.
    Blind to a one-row conservation defect of P_cell (arm C); see the
    module docstring."""
    result, params = _solve(cell)
    k_inf, _ = _infinite_medium(cell.groups)
    bound = 10 * params.keff_tol
    assert abs(result.keff - k_inf) <= bound, (
        f"{cell.id}: k = {result.keff!r}, k∞ = {k_inf!r}, "
        f"|Δk| = {abs(result.keff - k_inf):.3e} > {bound:.1e}"
    )


@pytest.mark.l1
@pytest.mark.verifies("cp-white-cell-infinite-medium", "p-inf", "cp-infinite-lattice-sum")
@pytest.mark.parametrize("cell", _SPATIAL_PARAMS)
def test_flux_is_flat_in_every_group(cell: WhiteCell) -> None:
    """C2: per group, ``max_i |φ_ig/⟨φ_g⟩ − 1| ≤ _FLATNESS_BAND``. The only
    row that sees a conservation defect of one row of P_cell (arm C), to
    which C1 and C3 are blind; itself blind to the scattering convention
    (arm D)."""
    result, _ = _solve(cell)
    deviation = np.max(np.abs(result.flux / _group_averages(cell, result) - 1.0), axis=0)
    assert np.all(deviation <= _FLATNESS_BAND), (
        f"{cell.id}: per-group flatness deviation {deviation} > {_FLATNESS_BAND:.0e}"
    )


@pytest.mark.l1
@pytest.mark.verifies("cp-white-cell-infinite-medium", "p-inf", "cp-infinite-lattice-sum")
@pytest.mark.parametrize("cell", _ALL_PARAMS)
def test_group_spectrum_is_the_infinite_medium_eigenvector(cell: WhiteCell) -> None:
    """C3: the normalised group spectrum equals the pencil's dominant
    eigenvector to ``10 × flux_tol``. Measured worst 4.9e-13
    (CYL-R0.5-uniform5-2g) against 1e-11. Sees the scattering-transpose
    convention (arm D); blind to arm C."""
    result, params = _solve(cell)
    _, reference = _infinite_medium(cell.groups)
    averages = _group_averages(cell, result)
    spectrum = averages / averages.sum()
    bound = 10 * params.flux_tol
    np.testing.assert_allclose(
        spectrum, reference, rtol=0.0, atol=bound,
        err_msg=f"{cell.id}: group spectrum is not the infinite-medium eigenvector",
    )


@pytest.mark.foundation
def test_the_509_cells_are_strict_xfails() -> None:
    """The #509 sphere cells carry ``xfail(strict=True)`` citing #509 in
    every parametrised row, and no other cell carries an xfail. A marker
    inside ``pytest.param(marks=…)`` that loses ``strict`` would turn a fix
    into a silent XPASS; this row makes that red. The #509 set is 3 cells
    (sphere R = 10 at 2G and 4G, sphere graded R = 4 at 4G)."""
    n_509 = 0
    for params in (_ALL_PARAMS, _SPATIAL_PARAMS):
        for param in params:
            (cell,) = param.values
            xfails = [m for m in param.marks if m.name == "xfail"]
            if _is_509(cell):
                n_509 += 1
                assert len(xfails) == 1, f"{cell.id}: expected one xfail, got {xfails}"
                assert xfails[0].kwargs.get("strict") is True, f"{cell.id}: xfail is not strict"
                assert "#509" in xfails[0].kwargs.get("reason", ""), f"{cell.id}: reason"
            else:
                assert not xfails, f"{cell.id}: unexpected xfail {xfails}"
    assert n_509 == 6, f"expected the 3 #509 cells in both lists (6 params), got {n_509}"


@pytest.mark.l1
def test_sphere_R10_flatness_defect_is_present() -> None:
    """RECORD (green today, red at the fix of #509): the 1G sphere,
    R = 10 cm, 10 uniform cells, with a white boundary is NOT flat; its
    flatness deviation is above 1e-3 (measured 4.8e-2 on 2026-09-25). At
    one group ``k = k∞`` exactly, so only the shape shows the defect. A
    partial fix that halves the deviation leaves this row green and the
    strict xfails failing; only a full fix turns this row red and the
    xfails into passes."""
    cell = WhiteCell(CoordSystem.SPHERICAL, 10.0, "uniform", 10, "1g")
    result, _ = _solve(cell)
    deviation = float(np.max(np.abs(result.flux / _group_averages(cell, result) - 1.0)))
    assert deviation > 1e-3, (
        f"#509 fixed: delete this row and the xfail marks (flatness deviation "
        f"now {deviation:.2e})"
    )
