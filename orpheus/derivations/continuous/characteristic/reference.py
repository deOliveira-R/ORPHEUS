r"""The characteristic reference: the Galerkin system posed from a specification, and the observables it answers.

:class:`CharacteristicDerivation` is the one door from the interface
vocabulary to the package. It reads a
:class:`~orpheus.specification.specification.GeometrySpecification` (the
materials, a finite geometry and one question) and a :class:`Resolution`,
refuses at construction what the reference does not serve, poses the
:class:`~.system.GalerkinSystem` the question needs, and answers every
observable of that question from one solve (the user's rulings of
2026-10-07, P1 step (b) fifth rung):

* :class:`~orpheus.numerics.question.Eigen` along the fission emission, at
  the physical point. Its fundamental mode is the k pencil's, its flux
  scaled so that the production the question declares
  (:attr:`~orpheus.numerics.question.Eigen.gauge`, by default fission and
  (n,2n), the functional the SN solver scales by) is 1 over the body. The mode nearest
  :math:`\tau` (:class:`~orpheus.numerics.question.Nearest`) is read in the
  k chart, the chart every reference reads
  :class:`~orpheus.numerics.observable.Eigenvalue` in (the charts of the
  parameters are owed to #529), and answers its eigenvalue only: a higher
  mode can produce no net fission neutrons (on a closed homogeneous body
  every one is biorthogonal to the flat adjoint), so the production gauge
  does not scale it, and nothing else is ruled;
* :class:`~orpheus.numerics.question.FixedSource`: the system posed with
  the source's regions, the least solution of its source pencil;
* :class:`~orpheus.numerics.question.Response`: the adjoint scalar flux
  :math:`R\psi^\dagger` of the detector, answered as the flux of the
  group-transposed problem
  (:meth:`~.cross_sections.RegionCrossSections.transposed`) with the
  detector as its source. For isotropic emission the transport is
  self-adjoint (reciprocity), so the adjoint flux of a problem is the
  forward flux of its adjoint cross sections, and every observable reads it
  with the forward machinery.

**The mesh-free functions enter by their role** (:mod:`orpheus.numerics.mesh_free_function`,
:mod:`orpheus.derivations.common.angular_measure`). The unknown is an
angle-integrated emission rate, so a source and a detector enter as the
retraction :math:`R` of their lift into phase space, and the :math:`4\pi`
is derived from the angular measure, never typed:

* a source table :math:`Q` lifts by the section :math:`E`, so its rate is
  :math:`R E Q = Q`;
* a detector table :math:`\Sigma_d` lifts by the pullback :math:`R^\dagger`,
  so the adjoint problem's rate is :math:`R R^\dagger \Sigma_d = 4\pi\Sigma_d`;
* a :class:`~orpheus.numerics.mesh_free_function.Symbolic` source or
  detector is already a density over the directions, so its rate is
  :math:`R f = \int f\,\mathrm d\Omega`;
* a weight pairs with the scalar flux as given.

Reciprocity then reads
:math:`\langle \Sigma_d, \phi(Q)\rangle = \langle Q, R\psi^\dagger(\Sigma_d)\rangle / 4\pi`
for two tables. A table is read onto the nodes exactly
(:meth:`~.basis.PanelBasis.on_nodes`); a symbolic function is projected
(:meth:`~.basis.PanelBasis.project`).

**A flux integral needs no point.** The pairing :math:`(W c_w)^{\mathsf T}
\phi_h`, with :math:`c_w` the coefficients of :math:`w` and :math:`\phi_h`
the Galerkin flux, is :math:`\int w\,\phi_h\,\mathrm dV`, since
:math:`\phi_h` is in the basis space and :math:`c_w` is the projection
:math:`P w`. The Galerkin flux satisfies :math:`W \phi_h = K q`, so it
reproduces every moment of the transported flux :math:`\mathcal K q`
against the basis, and the reading equals :math:`\langle w, \mathcal K
q\rangle` exactly for a weight in the basis space (every region-wise
constant); otherwise it differs from it by :math:`\langle w - P w,
\mathcal K q\rangle`, the weight's projection error against the
transported flux, which is the flux's projection error against the
weight. The flux at a point is the fifth rung's second half (5b).
"""

from __future__ import annotations

import operator
from collections.abc import Callable
from dataclasses import dataclass, field
from functools import cached_property
from typing import TYPE_CHECKING, NoReturn, assert_never

import numpy as np

from orpheus.data.cells import CellCoefficient, Channel
from orpheus.derivations.common.eigenvalue import production_emission
from orpheus.geometry.chord import ConcentricPartition
from orpheus.numerics.content import ContentIdentity
from orpheus.numerics.mesh_free_function import MeshFreeFunction, RegionwiseConstant, Symbolic
from orpheus.numerics.observable import Eigenvalue, FluxIntegral, Linear, PointValue
from orpheus.numerics.question import Eigen, FixedSource, Fundamental, Nearest, Response
from orpheus.numerics.traced_memo import traced_memo
from orpheus.reference.reading import Uncertified
from orpheus.reference.solution import ReferenceSolution
from orpheus.specification.specification import GeometrySpecification

from .assembly import TransportResolution
from .basis import PanelBasis
from .cross_sections import RegionCrossSections
from .system import GalerkinSystem
from .walls import Walls

if TYPE_CHECKING:
    import sympy

Lift = Callable[["sympy.Expr"], "sympy.Expr"]
"""An arrow from a per-region table into phase space: the section of a source, the pullback of a detector."""


def _section(rate: "sympy.Expr") -> "sympy.Expr":
    """A source's lift (:func:`~orpheus.derivations.common.angular_measure.section`)."""
    from orpheus.derivations.common.angular_measure import section

    return section(rate)


def _pullback(response: "sympy.Expr") -> "sympy.Expr":
    """A detector's lift (:func:`~orpheus.derivations.common.angular_measure.pullback`)."""
    from orpheus.derivations.common.angular_measure import pullback

    return pullback(response)


def _retraction(density: "sympy.Expr") -> "sympy.Expr":
    """The rate of a density over the directions (:func:`~orpheus.derivations.common.angular_measure.retraction`)."""
    from orpheus.derivations.common.angular_measure import retraction

    return retraction(density)


@dataclass(frozen=True, eq=False)
class Resolution(ContentIdentity):
    r"""The resolution of a characteristic reference: one value for the solve and every reading.

    Attributes
    ----------
    degree:
        The panel basis's polynomial degree :math:`p`.
    layers, ratio:
        The panels' geometric grading toward every wall and interface
        (:meth:`~.basis.PanelBasis.of`).
    transport:
        The transport blocks' line and traversal rules.
    source_points:
        The Gauss points per piece that project a symbolic source, detector
        or weight onto the basis (:meth:`~.basis.PanelBasis.project`); at
        least :math:`2(p + 1)`, so that a function of the panel space is
        projected onto itself.
    """

    degree: int
    layers: int
    ratio: float
    transport: TransportResolution
    source_points: int

    def __post_init__(self) -> None:
        for name in ("degree", "layers", "source_points"):
            try:
                object.__setattr__(self, name, operator.index(getattr(self, name)))
            except TypeError:
                raise TypeError(f"the {name} of a resolution is an integer, got {getattr(self, name)!r}") from None
        if not isinstance(self.transport, TransportResolution):
            raise TypeError(f"the transport resolution is a TransportResolution, got a {type(self.transport).__name__}")
        if self.source_points < 2 * (self.degree + 1):
            raise ValueError(
                f"a projection onto degree {self.degree} panels takes at least {2 * (self.degree + 1)} points per piece; "
                f"got {self.source_points}"
            )


# ── the answers ──────────────────────────────────────────────────────────


@dataclass(frozen=True)
class _FundamentalAnswer:
    r"""The fundamental mode: :math:`k` and the flux coefficients ``(G, N)``, total production 1."""

    k: float
    flux: np.ndarray


@dataclass(frozen=True)
class _ModeAnswer:
    r"""A mode nearest :math:`\tau`: its eigenvalue alone."""

    k: float


@dataclass(frozen=True)
class _SourceAnswer:
    """A fixed source's flux, or a detector's adjoint scalar flux, ``(G, N)``."""

    flux: np.ndarray


_Answer = _FundamentalAnswer | _ModeAnswer | _SourceAnswer


def _refuse_direction(function: MeshFreeFunction, role: str) -> None:
    r"""Refuse a symbolic source or detector that depends on the direction.

    SCOPE-BOUNDARY[guard] machinery: an anisotropic source or detector (its own first-flight transport along each line).
    ruling: the user, 2026-10-07, P1 step (b) fourth rung, Q2 (the emission is isotropic); the fifth rung's sketch.
    revisit: when a consumer poses an anisotropic source or detector to a reference.
    """
    if isinstance(function, Symbolic) and function.depends_on(Symbolic.mu, Symbolic.phi):
        raise NotImplementedError(
            f"the characteristic reference serves an isotropic {role} only: the {role} depends on the direction"
        )


def _refuse_point_reading() -> NoReturn:
    r"""Refuse the flux at a point.

    SCOPE-BOUNDARY[guard] machinery: the reading at a point, the transport of the converged emission over :meth:`~orpheus.geometry.chart.Chart.directions_at`.
    ruling: the user, 2026-10-07, P1 step (b) fifth rung, Q1: the rung splits, and the reading is 5b.
    revisit: P1 step (b), rung 5b.
    """
    raise NotImplementedError("the characteristic reference does not read the flux at a point yet (P1 step (b), rung 5b)")


def _refuse_mode_flux(k: float) -> NoReturn:
    r"""Refuse the flux of a mode nearest :math:`\tau`.

    SCOPE-BOUNDARY[guard] machinery: the scale of a higher mode (a gauge that does not vanish on it), owed to the eigen answer's contract.
    ruling: the user, 2026-10-07, P1 step (b) fifth rung, after the reviews: a Nearest answer reads its eigenvalue only.
    revisit: when a consumer reads a higher mode's flux, with #529's eigen-answer contract.
    """
    raise NotImplementedError(
        f"a higher mode has no flux scale (k = {k:.6g}): the production gauge can vanish on it, so a Nearest answer "
        "reads its eigenvalue only"
    )


@dataclass(frozen=True, eq=False)
class CharacteristicDerivation(ContentIdentity):
    """The characteristic reference of a specification at a resolution, solved on the first evaluation and read uncertified.

    Construction is the one door for what the posing decides. It refuses an
    infinite medium (no geometry); a question off the physical point; an
    :class:`~orpheus.numerics.question.Eigen` along any parameter but the
    fission emission; a source or a detector that depends on the direction;
    and, through the objects it builds, a mixture with anisotropic emission
    (:meth:`~.cross_sections.RegionCrossSections.of`), a wall law the
    reference does not read (:meth:`~.walls.Walls.of`) and a grading
    :meth:`~.basis.PanelBasis.of` refuses. What only the solve decides is
    refused at the first evaluation: a body supercritical for its fixed
    source, a nearest eigenvalue that is complex, a pencil with no
    fundamental mode. An observable the question does not answer (an
    eigenvalue of a source question, a point) is refused before any solve.

    Its content is its two fields; the system, the basis and the cross
    sections are derived from them, so :meth:`evaluate` is a traced memo
    keyed on the derivation and the observable (#405 P3).
    """

    specification: GeometrySpecification
    resolution: Resolution
    cross_sections: RegionCrossSections = field(init=False, repr=False, compare=False)
    basis: PanelBasis = field(init=False, repr=False, compare=False)
    system: GalerkinSystem = field(init=False, repr=False, compare=False)
    source: np.ndarray | None = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        specification, resolution = self.specification, self.resolution
        if not isinstance(specification, GeometrySpecification):
            raise TypeError(
                f"characteristic_reference answers a GeometrySpecification (lines need a finite body), "
                f"got a {type(specification).__name__}: an infinite medium has no geometry"
            )
        if not isinstance(resolution, Resolution):
            raise TypeError(f"the resolution is a Resolution, got a {type(resolution).__name__}")
        geometry, question = specification.geometry, specification.question
        if question.point:
            raise ValueError(f"characteristic_reference answers at the physical point, got the point {dict(question.point)!r}")
        walls = Walls.of(geometry)
        basis = PanelBasis.of(ConcentricPartition.of(geometry), resolution.degree, resolution.layers, resolution.ratio)
        object.__setattr__(self, "basis", basis)
        cross_sections = RegionCrossSections.of([specification.materials[m] for m in geometry.mat_ids])
        match question:
            case Eigen(parameter=parameter):
                fission = CellCoefficient.every(Channel.FISSION_EMISSION).resolve(specification.materials)
                if parameter != fission:
                    raise ValueError(
                        f"characteristic_reference answers an eigenvalue along the fission emission, got the parameter {parameter!r}"
                    )
                source = None
            case FixedSource(source=function):
                _refuse_direction(function, "source")
                source = self._coefficients(function, _section)
            case Response(detector=function):
                _refuse_direction(function, "detector")
                cross_sections = cross_sections.transposed()
                source = self._coefficients(function, _pullback)
            case unreachable:
                assert_never(unreachable)
        regions = None if source is None else self._regions_of(source)
        object.__setattr__(self, "cross_sections", cross_sections)
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "system", GalerkinSystem(basis, walls, cross_sections, resolution.transport, regions))

    # ── the functions on the basis ───────────────────────────────────────

    def _coefficients(self, function: MeshFreeFunction, lift: Lift | None) -> np.ndarray:
        r"""A mesh-free function's coefficients on the basis, ``(G, N)``.

        With a ``lift`` (a source's section, a detector's pullback), the rate
        the transport reads, the retraction of the function's lift into phase
        space; with none (a weight), the function as given.
        """
        match function:
            case RegionwiseConstant(values=values):
                scale = 1.0 if lift is None else float(_retraction(lift(_unit())))
                return scale * self.basis.on_nodes(values).T
            case Symbolic():
                return self._projection(function if lift is None else Symbolic.of(*map(_retraction, function.expressions)))
            case _:
                assert_never(function)

    def _projection(self, function: Symbolic) -> np.ndarray:
        """A symbolic function of the position alone, projected onto the basis, ``(G, N)``."""
        import sympy

        scalar = function.without(Symbolic.mu, Symbolic.phi)
        evaluators = [sympy.lambdify(Symbolic.r, expression, "numpy") for expression in scalar.expressions]

        def values(c: np.ndarray) -> np.ndarray:
            return np.stack([np.broadcast_to(np.asarray(f(c), dtype=float), c.shape) for f in evaluators])

        breakpoints = self.specification.geometry.breakpoints
        steps = scalar.steps((float(breakpoints[0]), float(breakpoints[-1])))
        return self.basis.project(values, self.resolution.source_points, steps)

    def _regions_of(self, coefficients: np.ndarray) -> np.ndarray:
        """The regions where each group's coefficients are not all zero, ``(n, G)``: where a source is posed."""
        regions = np.zeros((self.basis.regions.n_regions, coefficients.shape[0]), dtype=bool)
        np.logical_or.at(regions, self.basis.region, (coefficients != 0.0).T)
        return regions

    def _pairing(self, weight: np.ndarray, flux: np.ndarray) -> float:
        r""":math:`\sum_g (W c_{w,g})^{\mathsf T} \phi_g`: a weight's coefficients paired with flux coefficients, both ``(G, N)``."""
        return float(np.sum((weight @ self.basis.mass) * flux))

    # ── the answer ───────────────────────────────────────────────────────

    @cached_property
    def answer(self) -> _Answer:
        """The question's answer, solved once."""
        system = self.system
        match self.specification.question:
            case Eigen(mode=Fundamental(), gauge=CellCoefficient() as gauge):
                mode = system.pencil.fundamental()
                flux = system.flux(mode.vector)
                materials, mat_ids = self.specification.materials, self.specification.geometry.mat_ids
                emission = np.stack([production_emission(gauge, m, materials[m]) for m in mat_ids])  # (n, G): the declared production
                return _FundamentalAnswer(mode.k, flux / self._pairing(self.basis.on_nodes(emission).T, flux))
            case Eigen(mode=Nearest(tau=tau)):
                return _ModeAnswer(self._nearest(tau))
            case FixedSource() | Response():
                return _SourceAnswer(system.fixed_source(self.source))
            case unreachable:
                assert_never(unreachable)

    def _nearest(self, tau: float) -> float:
        r"""The eigenvalue of the k pencil nearest :math:`\tau`, among the modes that emit fission neutrons.

        A production matrix :math:`F K` of rank below its size (a group
        fission does not emit into, a region that does not produce) gives
        the pencil a cluster of eigenvalues at :math:`k = 0`, each a vector
        :math:`F K` annihilates, which rounding scatters around zero, some
        complex. They are no modes of the fission problem: a mode is kept
        when :math:`\lVert F K v\rVert` exceeds the rank tolerance of
        :math:`F K` (its size times the machine epsilon times its norm).
        """
        production = self.system.pencil.production
        tolerance = max(production.shape) * np.finfo(float).eps * np.linalg.norm(production, 2)
        modes = [m for m in self.system.pencil.spectrum() if np.linalg.norm(production @ m.vector) > tolerance]
        mode = min(modes, key=lambda m: abs(m.eigenvalue - tau))
        if mode.eigenvalue.imag != 0.0:
            raise ValueError(f"the eigenvalue nearest {tau:.6g} is complex, {mode.eigenvalue:.6g}")
        return float(mode.eigenvalue.real)

    @traced_memo
    def evaluate(self, observable: Eigenvalue | Linear) -> Uncertified:
        """The observable's value, uncertified: the family derives no bound yet (#566)."""
        match observable, self.specification.question:
            case Eigenvalue(), FixedSource() | Response():
                raise ValueError("an eigenvalue is read from an eigen question's answer")
            case PointValue(), _:
                _refuse_point_reading()
        match observable, self.answer:                                   # what remains needs the solve
            case Eigenvalue(), _FundamentalAnswer(k=k) | _ModeAnswer(k=k):
                return Uncertified(k)
            case FluxIntegral(weight=weight), _FundamentalAnswer(flux=flux) | _SourceAnswer(flux=flux):
                return Uncertified(self._pairing(self._coefficients(weight, None), flux))
            case FluxIntegral(), _ModeAnswer(k=k):
                _refuse_mode_flux(k)
        raise AssertionError(f"unreachable: {observable!r} on {type(self.answer).__name__}")


def _unit() -> "sympy.Expr":
    """The rate 1, which a lift and the retraction scale."""
    import sympy

    return sympy.Integer(1)


def characteristic_reference(specification: GeometrySpecification, resolution: Resolution) -> ReferenceSolution:
    r"""The characteristic reference's :class:`~orpheus.reference.solution.ReferenceSolution`, uncertified and lazy.

    The refusals of the posing are the derivation's
    (:class:`CharacteristicDerivation`), all before any solve: nothing is
    solved here.
    """
    return ReferenceSolution(specification, CharacteristicDerivation(specification, resolution), None)


__all__ = ["CharacteristicDerivation", "Resolution", "characteristic_reference"]
