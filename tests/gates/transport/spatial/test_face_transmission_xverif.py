r"""L1 cross-check: production's face transmission against the algebra of record.

Branch 1 is :mod:`orpheus.derivations.discrete.sn.face_transmission` (SymPy,
derived from each closure's balance and closure rule; gated by
``test_face_transmission_symbolic.py``).  Branch 2 is production numpy, reached
three ways, each a different code path:

1. **The cell kernel**: :func:`~orpheus.transport.spatial.scheme._face_transmission_matrix`
   drives a scheme's ``cell_kernel_batch`` one unit inflow at a time, at the
   streaming coefficients :math:`g_a = |\mu_a|/\Delta_a` (the kernel's
   ``s_axes``) and the total cross section :math:`\Sigma_t`.  Branch 1
   measures lengths in mean free paths, so production at
   :math:`(g, \Sigma_t)` is Branch 1 at :math:`g/\Sigma_t`.
1b. **The other blocks**: the same kernel's cell state and outflow under a
   unit inflow or a unit source moment, against ``inflow_to_cell``,
   ``source_to_cell`` and ``escape`` (no ``verifies``: the label states the
   transmission only).
2. **The verdict**: :meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.face_transmission_spectrum`
   (the gauge predicate's closure half, #344), which reads the kernel at two
   probe cells.
3. **The sweep**: one production sweep on a one-cell mesh, through
   :func:`orpheus.sn.loss_representation.default_for`, with a unit inflow on
   one face moment of one first-octant ordinate and the outflow read back from
   the :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`.
   `[M]` 2026-09-23 this is NOT the kernel path at every d: diamond at
   :math:`d = 1` runs ``CumprodScan`` (``affine_scan_coefficients``), at
   :math:`d = 2` ``ScanMarch`` (``cartesian_scan_coefficients``), and only at
   :math:`d = 3` the ``cell_kernel_batch`` wavefront; LD at :math:`d = 1`
   runs ``CumprodScan``, at :math:`d = 2` the kernel.  Legs 1 and 3 are
   therefore not redundant: the kernel and ``ScanMarch`` share the ``2 g``
   couplings of ``_cartesian_streaming_diagonal``, but ``ScanMarch`` spells
   its transmission ``a = 2 diag / S - 1`` itself and the 1-D scan spells its
   own ``2 |mu| A``, so a closure change on one path leaves the others alone
   (`[M]` the mutation arms P1a, P1b, P1c in the spec).

Independence (X4): the two branches share the balance equation they both
discretise (the bifurcation point) and nothing below it.  Diamond's Branch 1
is SymPy's ``linear_eq_to_matrix`` on the balance; production is a closed-form
numpy update.  LD's Branch 1 is ``ld_ubld`` (SymPy) and production is ``_ubld``
(numpy): two transcriptions of one Galerkin weak form, so LD's legs rest on the
independent grounds of that form (the Padé order of its 1-D transmission and
exactness on a bilinear flux), declared in ``rests_on``.

Every Branch-1 value is the module's own block evaluated at a point
(``cell_response(...).transmission(at=...)`` and its siblings, exact over QQ),
so these rows pin the module's composition as well as production.

Tolerance, per law: the transmission is one :math:`m \times m` cell solve
followed by a trace (:math:`m = 1` diamond, :math:`2^d` LD), so the backward
error bound is :math:`m\,\kappa(A)\,\varepsilon` relative to
:math:`\max(1, \lVert T\rVert_\infty)`; the gate is ten times that.  `[M]`
2026-09-23 over 6 points x 5 (closure, d) rows the kernel leg reads at most
1.0 of that unit (diamond :math:`d = 1`), and the sweep leg at most 1.0.

Markers: ``l1`` (an exact algebraic reference against production).  The
diamond rows state the label ``dd-face-transmission-spectrum`` and carry
``verifies``; LD's rows carry none (no label states them).
"""
from __future__ import annotations

import numpy as np
import pytest
import sympy as sp

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.derivations.discrete.sn.face_transmission import (
    LAMBDA,
    Closure,
    cell_response,
    streaming_coefficients,
    transmission_spectrum,
)
from orpheus.geometry import BC
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn import loss_representation
from orpheus.sn.problem import SNProblem
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.mesh.axis import AxisMesh
from orpheus.transport.spatial.diamond import DiamondDifference
from orpheus.transport.spatial.linear_discontinuous import LinearDiscontinuous
from orpheus.transport.spatial.scheme import (
    _PROBE_CELLS,
    FaceModeDamping,
    _face_transmission_matrix,
)

LABEL = "dd-face-transmission-spectrum"
_SYM = "tests/gates/transport/spatial/test_face_transmission_symbolic.py"
_UBLD = "tests/gates/transport/spatial/test_ld_ubld_symbolic.py"
_EPS = float(np.finfo(float).eps)
_SAFETY = 10.0

_SCHEME = {
    Closure.DIAMOND: DiamondDifference,
    Closure.LINEAR_DISCONTINUOUS: LinearDiscontinuous,
}

R = sp.Rational
#: (g per axis, Sigma_t): generic; symmetric thin; strongly anisotropic;
#: symmetric thick; near-void (tau ~ 1e-6); extreme thick (tau ~ 1e6).
_POINTS = (
    ((R(7, 10), R(13, 10), R(2, 5)), R(9, 10)),
    ((R(1), R(1), R(1)), R(3, 10)),
    ((R(1, 100), R(50), R(1)), R(1)),
    ((R(1), R(1), R(1)), R(50)),
    ((R(1), R(2), R(3)), R(1, 10**6)),
    ((R(1, 1000), R(1, 1000), R(1, 1000)), R(1000)),
)
_POINT_IDS = ("generic", "thin", "anisotropic", "thick", "near-void", "extreme-thick")


def _fail_unless(condition: bool, message: str) -> None:
    """``pytest.fail`` unless ``condition`` holds; survives ``python -O``."""
    if not condition:
        pytest.fail(message)


def _branch1_at(closure: Closure, d: int, g) -> tuple[np.ndarray, float]:
    r"""Branch 1's :math:`T` at the point ``g`` (the module's ``transmission(at=)``, exact over QQ), and :math:`\kappa(A)` there."""
    r = cell_response(closure, d)
    at = dict(zip(r.g, g))
    return (
        np.array(r.transmission(at).evalf(30), dtype=float),
        float(np.linalg.cond(np.array(r.cell_operator.subs(at).evalf(30), dtype=float))),
    )


def _assert_transmission(produced: np.ndarray, closure: Closure, d: int, g, where: str) -> None:
    expected, kappa = _branch1_at(closure, d, g)
    _fail_unless(
        produced.shape == expected.shape,
        f"{where}: shape {produced.shape}, Branch 1 {expected.shape} "
        f"(a face reduced to its first moment reads as a plausible smaller map)",
    )
    m = cell_response(closure, d).n_cell_moments
    bound = _SAFETY * m * kappa * _EPS * max(1.0, float(np.abs(expected).max()))
    gap = float(np.abs(produced - expected).max())
    _fail_unless(
        gap <= bound,
        f"{where}: max|T_production - T_branch1| = {gap:.3e} > {bound:.3e} "
        f"(= {_SAFETY:g} x {m} x kappa {kappa:.2f} x eps); max|T| = {np.abs(expected).max():.3e}",
    )


_KERNEL_ROWS = [
    pytest.param(Closure.DIAMOND, d, marks=pytest.mark.verifies(LABEL), id=f"diamond-{d}")
    for d in (1, 2, 3)
] + [pytest.param(Closure.LINEAR_DISCONTINUOUS, d, id=f"ld-{d}") for d in (1, 2)]


# ── leg 1: the cell kernel ───────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.verifies(LABEL)
@pytest.mark.rests_on(
    f"{_SYM}::test_diamond_transmission_is_the_page_closed_form",
    f"{_SYM}::test_stability_function_is_an_a_stable_pade_approximant",
)
@pytest.mark.parametrize(("g", "sigma_t"), _POINTS, ids=_POINT_IDS)
@pytest.mark.parametrize("d", (1, 2, 3))
def test_diamond_kernel_applies_the_branch1_transmission(d, g, sigma_t) -> None:
    r"""``DiamondDifference.cell_kernel_batch`` induces :math:`(2/D)\mathbf 1 w^{\mathsf T} - I`."""
    produced = _face_transmission_matrix(
        DiamondDifference(), d, tuple(float(x) for x in g[:d]), float(sigma_t),
    )
    _assert_transmission(produced, Closure.DIAMOND, d, [x / sigma_t for x in g[:d]], f"diamond kernel d={d}")


@pytest.mark.l1
@pytest.mark.rests_on(
    f"{_SYM}::test_stability_function_is_an_a_stable_pade_approximant",
    f"{_SYM}::test_step_and_ld_damp_every_face_mode",
    f"{_UBLD}::TestOracleId1Reduction::test_d1_reduction_to_production_schur",
    f"{_UBLD}::TestOracleIIBilinearExactness::test_d2_exact_on_bilinear",
)
@pytest.mark.parametrize(("g", "sigma_t"), _POINTS, ids=_POINT_IDS)
@pytest.mark.parametrize("d", (1, 2))
def test_ld_kernel_applies_the_branch1_transmission(d, g, sigma_t) -> None:
    r"""``LinearDiscontinuous.cell_kernel_batch`` induces Branch 1's :math:`d\,2^{d-1}`-square map.

    :math:`d = 3` is absent because production cannot drive it: its numpy
    inflow lift handles ``axis in {0, d-1}`` only (see the strict xfail below).
    """
    produced = _face_transmission_matrix(
        LinearDiscontinuous(), d, tuple(float(x) for x in g[:d]), float(sigma_t),
    )
    _assert_transmission(
        produced, Closure.LINEAR_DISCONTINUOUS, d, [x / sigma_t for x in g[:d]], f"LD kernel d={d}",
    )


# ── leg 1b: the other three blocks of the cell response ──────────────────


def _drive_kernel(scheme, d: int, g, sigma_t, *, inflow=None, source=None):
    """One ``cell_kernel_batch`` call: a unit inflow moment ``(axis, k)`` or a unit source moment ``k``.

    Returns the cell moments and the outflow face moments, flattened.
    """
    per_face = scheme.spatial_basis_per_axis ** (d - 1)
    face_shape = (1, 1, 1, per_face) if per_face > 1 else (1, 1, 1)
    psi_in = []
    for a in range(d):
        face = np.zeros(face_shape)
        if inflow is not None and inflow[0] == a:
            face.reshape(-1)[inflow[1]] = 1.0
        psi_in.append(face)
    moments = scheme.spatial_basis_per_axis**d
    q = np.zeros((1, 1, 1, moments) if moments > 1 else (1, 1, 1))
    if source is not None:
        q.reshape(-1)[source] = 1.0
    cell, out = scheme.cell_kernel_batch(
        psi_in=tuple(psi_in), s_axes=tuple(np.full((1, 1, 1), float(x)) for x in g),
        reaction_xs=np.full((1, 1), float(sigma_t)), Q_cells=q,
    )
    return np.ravel(cell), np.concatenate([np.ravel(o) for o in out])


@pytest.mark.l1
@pytest.mark.rests_on(
    f"{_SYM}::test_particle_balance_from_the_four_blocks",
    f"{_SYM}::test_flat_flux_passes_through_the_cell",
    "tests/gates/transport/spatial/test_face_transmission_xverif.py::test_diamond_kernel_applies_the_branch1_transmission",
    "tests/gates/transport/spatial/test_face_transmission_xverif.py::test_ld_kernel_applies_the_branch1_transmission",
)
@pytest.mark.parametrize(("g", "sigma_t"), _POINTS, ids=_POINT_IDS)
@pytest.mark.parametrize(
    ("closure", "d"),
    [(Closure.DIAMOND, d) for d in (1, 2, 3)] + [(Closure.LINEAR_DISCONTINUOUS, d) for d in (1, 2)],
    ids=["diamond-1", "diamond-2", "diamond-3", "ld-1", "ld-2"],
)
def test_kernel_applies_the_branch1_escape_and_cell_blocks(closure, d, g, sigma_t) -> None:
    r"""The kernel's cell state and outflow equal ``inflow_to_cell``, ``source_to_cell`` and ``escape``.

    Unit conventions `[M]` 2026-09-23: production's source ``Q`` is per unit
    volume in units of :math:`\Sigma_t\psi`, so Branch 1's source blocks
    (at :math:`\Sigma_t = 1`) are divided by :math:`\Sigma_t`; the inflow
    block is dimensionless.  For diamond the cell state is the average; for LD
    it is the :math:`2^d` Kronecker moment vector, source moments in the same order.
    """
    scheme = _SCHEME[closure]()
    r = cell_response(closure, d)
    at = dict(zip(r.g, [x / sigma_t for x in g[:d]]))
    gd = g[:d]
    blocks = {
        "inflow_to_cell": np.array(r.inflow_to_cell(at).evalf(30), dtype=float),
        "source_to_cell": np.array(r.source_to_cell(at).evalf(30), dtype=float) / float(sigma_t),
        "escape": np.array(r.escape(at).evalf(30), dtype=float) / float(sigma_t),
    }
    per_face = r.n_face_moments // d
    to_cell = np.array([
        _drive_kernel(scheme, d, gd, sigma_t, inflow=divmod(j, per_face))[0] for j in range(r.n_face_moments)
    ]).T
    source_runs = [_drive_kernel(scheme, d, gd, sigma_t, source=k) for k in range(r.source_mass.cols)]
    produced = {
        "inflow_to_cell": to_cell,
        "source_to_cell": np.array([cell for cell, _ in source_runs]).T,
        "escape": np.array([out for _, out in source_runs]).T,
    }
    kappa = float(np.linalg.cond(np.array(r.cell_operator.subs(at).evalf(30), dtype=float)))
    for name, expected in blocks.items():
        got = produced[name]
        _fail_unless(got.shape == expected.shape, f"{closure.value} d={d} {name}: shape {got.shape} vs {expected.shape}")
        bound = _SAFETY * r.n_cell_moments * kappa * _EPS * max(1.0, float(np.abs(expected).max()))
        gap = float(np.abs(got - expected).max())
        _fail_unless(gap <= bound, f"{closure.value} d={d} {name}: gap {gap:.3e} > {bound:.3e}")


# ── leg 2: the verdict and the spectral radius ───────────────────────────


def _branch1_radius(closure: Closure, d: int) -> float:
    """Max |root| of Branch 1's factors over production's two probe cells."""
    spectrum = transmission_spectrum(closure, d)
    radii = []
    for streaming, sigma_t in _PROBE_CELLS:
        at = {
            ga: sp.nsimplify(streaming[a]) / sp.nsimplify(sigma_t)
            for a, ga in enumerate(streaming_coefficients(d))
        }
        roots = np.concatenate([
            np.roots([complex(c) for c in sp.Poly(fac.as_expr().subs(list(at.items())), LAMBDA).all_coeffs()])
            for fac, _ in spectrum.factors
        ])
        radii.append(float(np.abs(roots).max()))
    return max(radii)


def _theorem_verdict(closure: Closure, d: int) -> FaceModeDamping:
    return (
        FaceModeDamping.DAMPED
        if transmission_spectrum(closure, d).damps_every_face_mode
        else FaceModeDamping.UNDAMPED
    )


def _assert_verdict(closure: Closure, d: int) -> None:
    produced = _SCHEME[closure]().face_transmission_spectrum(d)
    _fail_unless(
        produced.damping is _theorem_verdict(closure, d),
        f"{closure.value}, d={d}: production says {produced.damping}, the theorem "
        f"{_theorem_verdict(closure, d)} ({produced.undetermined_because})",
    )
    radius = _branch1_radius(closure, d)
    _fail_unless(
        produced.spectral_radius is not None
        and abs(produced.spectral_radius - radius) <= 16 * _EPS,
        f"{closure.value}, d={d}: production rho = {produced.spectral_radius!r}, "
        f"Branch 1 at the probe cells {radius!r}",
    )


@pytest.mark.l1
@pytest.mark.rests_on(
    f"{_SYM}::test_diamond_carries_exactly_d_minus_1_undamped_face_modes",
    f"{_SYM}::test_step_and_ld_damp_every_face_mode",
    "tests/gates/transport/spatial/test_face_transmission_xverif.py::test_diamond_kernel_applies_the_branch1_transmission",
    "tests/gates/transport/spatial/test_face_transmission_xverif.py::test_ld_kernel_applies_the_branch1_transmission",
)
@pytest.mark.parametrize(("closure", "d"), _KERNEL_ROWS)
def test_face_transmission_spectrum_agrees_with_the_theorem(closure, d) -> None:
    r"""Production's DAMPED / UNDAMPED verdict and :math:`\rho` equal Branch 1's.

    `[M]` 2026-09-23 the radii agree to at most 1 ULP (diamond reads
    ``1.0000000000000002``); the tolerance is 16 ULP.  A spectrum is blind to
    a transpose or a similarity of :math:`T` (Mode 12), which is why the
    kernel leg compares the matrix and this one rests on it.
    """
    _assert_verdict(closure, d)


@pytest.mark.l1
@pytest.mark.rests_on(f"{_SYM}::test_step_and_ld_damp_every_face_mode")
@pytest.mark.xfail(
    strict=True,
    reason=(
        "#503: production LD cannot be driven at ndim=3: orpheus.transport.spatial._ubld."
        "assemble_inflow_axis handles axis in {0, d-1} only, so face_transmission_spectrum "
        "reads UNDETERMINED; Branch 1 proves LD at d=3 DAMPED, rho = 0.9048361934477379 at "
        "the probe cells. The paired record is test_face_transmission_damping.py::"
        "test_UNDETERMINED_is_a_third_state_and_it_carries_its_reason"
    ),
)
def test_ld_d3_verdict_agrees_with_the_theorem() -> None:
    r"""The challenge row: production LD at :math:`d = 3` must read DAMPED, :math:`\rho` = Branch 1's.

    Red today for the named reason (``UNDETERMINED``); it turns green, and so
    fails as a strict XPASS, the day production gains the interior-axis lift,
    together with its paired record row, which then reddens.
    """
    _assert_verdict(Closure.LINEAR_DISCONTINUOUS, 3)


@pytest.mark.foundation
def test_the_ld_d3_challenge_row_is_strict() -> None:
    """The xfail above is strict; a mark moved into ``pytest.param`` would lose it."""
    marks = [m for m in test_ld_d3_verdict_agrees_with_the_theorem.__dict__["pytestmark"] if m.name == "xfail"]
    _fail_unless(len(marks) == 1 and marks[0].kwargs.get("strict") is True, f"xfail marks: {marks}")


# ── leg 3: the production sweep on a one-cell mesh ───────────────────────

#: Two pure absorbers, so one sweep carries two different optical thicknesses
#: and a group-index slip (Mode 2) moves the answer.
_SIGMA_T = np.array([0.9, 3.7])
_WIDTHS = (0.7, 1.3, 0.4)
_AXES = "xyz"


def _one_cell_sweep_maps(scheme, d: int):
    """Per group: the octant map the production sweep applies, and its g."""
    mixture = make_mixture(
        sig_t=_SIGMA_T, sig_c=_SIGMA_T.copy(), sig_f=np.zeros(2), nu=np.zeros(2),
        chi=np.zeros(2), sig_s=np.zeros((2, 2)),
    )
    vacuum = BC("vacuum")
    axes = tuple(
        AxisMesh(edges=np.array([0.0, w]), bc_low=vacuum, bc_high=vacuum) for w in _WIDTHS[:d]
    )
    quadrature = Quadrature.level_symmetric(4) if d > 1 else Quadrature.gauss_legendre(n_ordinates=4)
    problem = SNProblem.from_axes(
        axes, quadrature, {0: mixture}, mat_map=np.zeros((1,) * d, dtype=int), scheme=scheme,
    )
    representation = loss_representation.default_for(problem, problem.scheme, problem.angular_closure)
    q = problem.quad
    cosines = (q.mu_x, q.mu_y, q.mu_z)[:d]
    first_octant = np.logical_and.reduce([mu > 0 for mu in cosines])
    n = int(np.argmax(first_octant))
    _fail_unless(bool(first_octant[n]), f"d={d}: no first-octant ordinate")
    stratum = representation.bind_sigma(problem.sigma_t_cell)
    cell_moments = (2**d,) if isinstance(scheme, LinearDiscontinuous) else ()
    ng = _SIGMA_T.size
    per_face = AngularBoundaryFlux.zeros(problem.angular_trace).face_view("xmin")[n].size // ng
    maps = []
    for group in range(ng):
        columns = []
        for a in range(d):
            for moment in range(per_face):
                flux = AngularBoundaryFlux.zeros(problem.angular_trace)
                flux.face_view(f"{_AXES[a]}min")[n].reshape(ng, per_face)[group, moment] = 1.0
                representation.sweep(np.zeros((q.N, ng) + (1,) * d + cell_moments), stratum, flux)
                columns.append(np.concatenate([
                    flux.face_view(f"{_AXES[b]}max")[n].reshape(ng, per_face)[group] for b in range(d)
                ]))
        g = [sp.Float(abs(float(cosines[a][n])) / (_WIDTHS[a] * _SIGMA_T[group]), 40) for a in range(d)]
        maps.append((np.array(columns).T, g))
    return maps, type(representation).__name__


_SWEEP_ROWS = [
    pytest.param(Closure.DIAMOND, 1, "CumprodScan", marks=pytest.mark.verifies(LABEL), id="diamond-1"),
    pytest.param(Closure.DIAMOND, 2, "ScanMarch", marks=pytest.mark.verifies(LABEL), id="diamond-2"),
    pytest.param(Closure.DIAMOND, 3, "FullFieldWavefront", marks=pytest.mark.verifies(LABEL), id="diamond-3"),
    pytest.param(Closure.LINEAR_DISCONTINUOUS, 1, "CumprodScan", id="ld-1"),
    pytest.param(Closure.LINEAR_DISCONTINUOUS, 2, "MovingFrontierWindow", id="ld-2"),
]


@pytest.mark.l1
@pytest.mark.rests_on(
    "tests/gates/transport/spatial/test_face_transmission_xverif.py::test_diamond_kernel_applies_the_branch1_transmission",
    "tests/gates/transport/spatial/test_face_transmission_xverif.py::test_ld_kernel_applies_the_branch1_transmission",
    f"{_SYM}::test_diamond_transmission_is_the_page_closed_form",
)
@pytest.mark.parametrize(("closure", "d", "representation"), _SWEEP_ROWS)
def test_production_sweep_applies_the_branch1_transmission(closure, d, representation) -> None:
    r"""One sweep through a one-cell mesh applies Branch 1's octant map, per group.

    The representation is asserted because the row's coverage claim is about
    that path: if ``default_for`` changes its choice, this docstring and the
    module docstring's path table change with it.
    """
    maps, ran = _one_cell_sweep_maps(_SCHEME[closure](), d)
    _fail_unless(ran == representation, f"d={d}: default_for chose {ran}, the row covers {representation}")
    for group, (produced, g) in enumerate(maps):
        _assert_transmission(produced, closure, d, g, f"{closure.value} sweep d={d} group {group}")


@pytest.mark.l1
@pytest.mark.xfail(
    strict=True,
    raises=NotImplementedError,
    reason="#503: production LD sweeps at d=3 raise in _ubld.assemble_inflow_axis (interior axis)",
)
def test_production_ld_sweep_at_d3_applies_the_branch1_transmission() -> None:
    """The d = 3 LD sweep, red today for the named reason only (``raises=``)."""
    maps, _ = _one_cell_sweep_maps(LinearDiscontinuous(), 3)
    for group, (produced, g) in enumerate(maps):
        _assert_transmission(produced, Closure.LINEAR_DISCONTINUOUS, 3, g, f"LD sweep d=3 group {group}")
