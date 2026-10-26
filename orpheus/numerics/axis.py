r"""Space-factor axes — the generators of axis-composed function spaces.

An **axis** is one tensor factor of a function space: the value object
recording what fixes that factor's coordinate space — its index shape, its
measure, its kind, and the manifold and atoms it was sampled from — together
with its **origin**, the object that minted it, which is provenance and
deliberately outside that identity.
Spaces are ordered products of axes (``FunctionSpace.of_axes``); the axis is
the unit the composition machinery reasons about — partitions, collapses,
frames, and (later) ⊕-lifts act **per axis**, never on an anonymous
position of a monolithic shape tuple.

The slots, precisely
====================

* ``label`` — the factor's role (``"angular"``, ``"spatial"``,
  ``"energy"``, …). Two factors that agree on everything else and play
  different roles are different axes.
* ``shape`` — the index set of this factor, rank ≥ 1. Rank > 1 is
  admissible (a spherical-harmonic axis is ``(L+1, 2L+1)``; a rank-``d``
  spatial axis is a legal design choice for CS2).
* ``weights`` — the **factor measure** over exactly ``shape``.
  ``None`` **is the counting measure, deliberately and always**: an axis
  has no "unbound" state, so the two-state ambiguity of the legacy
  ``FunctionSpace.inner_product_weights`` (``None`` = "no canonical
  quadrature" *or* "Euclidean") cannot arise on this type.
* ``origin`` — **provenance, never identity**: a closed sum of five arms
  (ruling 64 of the numerics untangling), read by exhaustive ``match``:

  ============================  =====================  ==================  =================
  arm                           kind it implies        ``support``         ``nodes``
  ============================  =====================  ==================  =================
  :class:`FromMeasure` (μ)      ``NODAL``              ``μ.support``       ``μ.nodes``
  :class:`FromQuadrature` (q)   ``NODAL``              ``q.measure``'s     ``q.measure``'s
  :class:`FromBasis` (b)        its layout's kind      ``b.domain``        none
  :class:`FromFrame` (f)        ``Modal(Φᵀw)``         ``f.basis.domain``  none
  :class:`Declared` (k)         ``k``                  none                none
  ============================  =====================  ==================  =================

  The generating object rides inside its arm so a consumer can recover
  what the axis forgot (direction cosines, the level fibration, the frame
  that can re-dress a head at another order) by matching the arm, and a
  consumer that needs a quadrature names that need in its own refusal.
  ``Declared`` is the honest arm where no generating object exists or may
  be imported here: the energy axes (their content is their edges), the
  homogeneous one-point axis and the rank-``d`` spatial arm (contract G6b;
  ruling 82). The LD scheme's spatial-moment axis is ``FromBasis``, the
  continuum head of its cell basis (#423).
* ``kind``, ``support``, ``nodes``, ``layout`` — **derived from the origin**, never
  passed: :attr:`Axis.kind` is the closed sum ``BasisKind.NODAL |``
  :class:`Modal` ``(unit)`` (ruling 41), so ``Axis(kind=NODAL,
  generator=basis)``, a nodal axis claiming a spectral generator, is
  unspellable. ``NODAL`` factors carry a coordinate cone (per-component
  positivity is meaningful); ``MODAL`` factors do not (a spectral
  coefficient may be negative for a positive function), and the modal arm
  CARRIES the factor's unit. ``support`` is the manifold the factor's
  coordinates live on, ``nodes`` the atoms a nodal factor samples (stored
  read-only), and ``layout`` how a basis indexes its coefficients (ruling
  81, :mod:`orpheus.numerics.coefficient_layout`): ``None`` on a sampled or
  declared axis, the basis's own layout on a ``FromBasis`` / ``FromFrame``
  head. The one generic coefficient axis is this class; the family is its
  layout and a Legendre head's spent axis is its domain (until S7 two
  subclasses, ``HarmonicAxis`` and ``LegendreAxis``, carried them).
* :attr:`Axis.unit` (read from the kind) — the coordinates of the
  constant function :math:`1` in this factor's own representation (ruling
  37: minted with the axis, by whatever mints it, like the weights). On a
  ``NODAL`` factor the coordinates are point values, so the unit is all
  ones. A modal factor is minted as ``Modal(unit)``: the coefficients of
  :math:`1` in the basis on a continuum head (:math:`e_0`, since
  :math:`Y_{00} = P_0 = 1`), the moments of :math:`1`,
  :math:`\Phi^{\mathsf T} w`, on a head a frame has dressed (its
  coordinates are covariant moments), :math:`e_0` on the LD scheme's
  spatial-moment axis. The marginal over the axis integrates against the
  covector :math:`m = G_k\,u` (the axis's metric block applied to the
  unit; :mod:`orpheus.numerics.axis_marginal`).

Canonical storage — one spelling per measure (ruled 2026-08-20)
===============================================================

Two construction rulings make the measure's identity unambiguous:

* **All-ones weights collapse to ``None`` at construction.** The counting
  measure has ONE spelling and therefore one identity; without this,
  ``weights=None`` and ``weights=np.ones(shape)`` would be the same
  measure with unequal identities — a twin exactly of the shape the
  fresh-``EnergyGrid``-per-access trap takes at the grid layer.
* **Weights are canonicalized as ``w + 0.0`` and stored read-only.**
  ``-0.0`` and ``+0.0`` are one measure and must be one identity/byte
  pattern; the addition also guarantees a defensive copy, so mutating the
  caller's array can never move an axis's hash after it has been used as
  a dict key. The atoms are stored the same way.

Non-finite weights are **refused** (a measure's weights are finite
numbers). There is deliberately **no non-negativity guard**: CS2's
quadrature axes legally carry signed weights (e.g. level-symmetric
families with negative weights), and the axis is the wrong layer to
outlaw them.

Identity — what fixes the factor's space (ruling 71)
====================================================

Equality and hash compare the content ``(type, label, shape, weights,
kind, support, nodes, layout)`` plus each subclass's own fields, never object
identity and never the origin. A nodal axis is therefore fixed by its
support, its atoms, its measure and its label; a modal one by its domain,
its kind (with the unit), its measure and its label. Two consequences
#631 measured:

* two rules with equal weights and different directions
  (``level_symmetric(4)`` against ``product(2, 12)``) are two axes, and the
  same μ-rule about two poles (``gauss_legendre_on_polar_orbit(8, "x")``
  against ``(8, "z")``) is two axes, separated by the support alone;
* the origin is provenance: ``Axis.of_quadrature(q)`` and
  ``Axis.of_measure(q.measure, "angular")`` are one axis.

Two axes that differ only in measure are **different axes** (the collapse
doctrine's "a genuine one-cell slab keeps its axis with weight V ≠ 1,
distinguished from the quotient point by measure"); an
:class:`EnergyAxis` never equals a generic :class:`Axis` carrying the same
field tuple. This identity is what ``FunctionSpace.of_axes`` derives
*space names* from, and what a space's identity key reads on an
axis-built space (``FunctionSpace.identity_key``).

Layering note: this module is model-independent numerics and names
nothing in ``orpheus.data``, not even for typing (the layer gate forbids
``numerics → data``, guarded or not). The energy axis is minted from its
CONTENT, its boundary energies (:meth:`EnergyAxis.of_edges`); the views
that read a nuclear-data object as an axis live above, in data: the grid's
axis view is ``EnergyGrid.as_axis`` and the energy-arm rule over a set of
mixtures is ``orpheus.data.materials.energy_axis_of``. The generating objects
(measure, quadrature, basis, frame) rank below the axis in the numerics
concept order and are read only through their attributes.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from enum import Enum, unique
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from .coefficient_layout import CellLayout, CoefficientLayout, HarmonicLayout, LegendreLayout, ModeLayout, PolynomialLayout, TensorLayout
from orpheus.values.content import ContentIdentity

if TYPE_CHECKING:
    from orpheus.numerics.basis.base import Basis, TruncatedBasis
    from orpheus.numerics.frame import FrameBase
    from orpheus.geometry.manifold import Manifold
    from orpheus.numerics.measure import DiscreteMeasure
    from orpheus.numerics.quadrature.directional import Quadrature

__all__ = [
    "Axis",
    "AxisKind",
    "BasisKind",
    "Declared",
    "EnergyAxis",
    "FromBasis",
    "FromFrame",
    "FromMeasure",
    "FromQuadrature",
    "Modal",
    "Origin",
]



@unique
class BasisKind(Enum):
    """The NODAL arm of a space factor's kind (the other arm is :class:`Modal`).

    ``NODAL`` — components are point/cell VALUES (indicator-like basis):
    the factor carries a coordinate cone, so per-component positivity is a
    meaningful predicate, and the coordinates of the constant function are
    all ones, so the arm needs no data. The MODAL arm carries data, the
    factor's unit, and is therefore the value class :class:`Modal`, not an
    enum member (ruling 41: the kind is the closed sum ``NODAL | Modal(unit)``).
    """

    NODAL = "nodal"

    def unit_over(self, shape: tuple[int, ...]) -> NDArray:
        r"""All ones over ``shape``: on a nodal factor the coordinates are
        point values, and the constant function's are ones."""
        ones = np.ones(shape)
        ones.setflags(write=False)
        return ones


@dataclass(frozen=True, eq=False)
class Modal(ContentIdentity):
    r"""The MODAL arm of a space factor's kind, carrying the factor's unit.

    Components are expansion COEFFICIENTS (a spectral basis): no coordinate
    cone, since a positive function may have negative coefficients, so
    per-component sign tests are meaningless and must be refused, not
    answered. A modal factor is minted WITH ``unit``, the coordinates of the
    constant function :math:`1` in its representation (rulings 37, 40 and
    41), so a modal factor that cannot integrate is unrepresentable. The
    unit is stored canonicalized (``u + 0.0``), copied and read-only, and
    must be finite; its shape is checked against the axis's shape by the
    axis. A ZERO unit is legal (the user's ruling, 2026-10-21, at S7; it
    amends ruling 41's "nonzero"): it is the honest unit of a factor whose
    span is orthogonal to the constants (the loss-kernel gauge's kernel
    modes), and the marginal over that factor integrates every field to 0,
    which is what each of its functions integrates to. A family whose span
    contains :math:`1` mints its unit nonzero by construction (:math:`e_0`,
    or the frame's moments of 1).
    """

    unit: NDArray

    def __post_init__(self) -> None:
        u = np.ascontiguousarray(self.unit, dtype=float)
        if not bool(np.isfinite(u).all()):
            raise ValueError(
                f"a modal unit must be finite (the coordinates of the constant "
                f"function 1); got {u!r}"
            )
        u = u + 0.0
        u.setflags(write=False)
        object.__setattr__(self, "unit", u)

    def unit_over(self, shape: tuple[int, ...]) -> NDArray:
        r"""The carried unit, which must live over exactly ``shape``."""
        if self.unit.shape != shape:
            raise ValueError(
                f"the modal unit must live over exactly the axis shape: unit "
                f"shape {self.unit.shape} != axis shape {shape}"
            )
        return self.unit


#: A space factor's kind, the closed sum of its two arms (ruling 41).
AxisKind = BasisKind | Modal




# ---------------------------------------------------------------------------
# The origin — what minted an axis (ruling 64)
# ---------------------------------------------------------------------------
#
# Each arm answers the three reads an axis derives its content from: the kind,
# the manifold the coordinates live on, and the atoms a nodal factor samples.
# An arm holding a generating object compares that object by IDENTITY
# (``eq=False``): the arm is provenance, and the content it answers is
# compared on the axis, where it is copied and read-only. ``Declared`` holds
# only a kind, so it is a value.


def _read_only(array: NDArray) -> NDArray:
    out = np.asarray(array, dtype=float) + 0.0
    out.setflags(write=False)
    return out


@dataclass(frozen=True)
class Declared:
    r"""The origin of an axis no generating object minted: the kind is declared.

    Used where the content is complete without atoms — the energy axes
    (their edges are their content), the homogeneous one-point axis and the
    rank-``d`` spatial arm (contract G6b: no rank-``d`` measure-to-axis
    pairing exists yet; ruling 82).
    """

    declared: AxisKind

    def __post_init__(self) -> None:
        if not isinstance(self.declared, (BasisKind, Modal)):
            raise TypeError(
                f"Declared: the kind must be BasisKind.NODAL or Modal(unit), "
                f"got {self.declared!r}"
            )

    @property
    def kind(self) -> AxisKind:
        return self.declared

    @property
    def support(self) -> None:
        return None

    @property
    def nodes(self) -> None:
        return None

    @property
    def layout(self) -> None:
        return None


@dataclass(frozen=True, eq=False)
class FromMeasure:
    r"""The origin of the axis a discrete measure samples: ``NODAL``, on the
    measure's support, at its atoms."""

    measure: DiscreteMeasure

    @property
    def kind(self) -> AxisKind:
        return BasisKind.NODAL

    @property
    def support(self) -> Manifold:
        return self.measure.support

    @property
    def nodes(self) -> NDArray:
        return self.measure.nodes

    @property
    def layout(self) -> None:
        return None


@dataclass(frozen=True, eq=False)
class FromQuadrature:
    r"""The origin of the angular axis a directional quadrature samples.

    The content is the rule's measure's (``NODAL``, its support, its
    atoms), so ``FromQuadrature(q)`` and ``FromMeasure(q.measure)`` mint one
    axis. The arm exists because the measure alone does not answer what the
    angular consumers recover through the space: the direction cosines and
    the level fibration (``level_indices``) live on the quadrature.
    ⚠ A LIVE REFERENCE: ``Quadrature`` is a mutable dataclass; the axis's
    own copies of the weights and atoms are the authoritative content.
    """

    quadrature: Quadrature

    @property
    def kind(self) -> AxisKind:
        return BasisKind.NODAL

    @property
    def support(self) -> Manifold:
        return self.quadrature.measure.support

    @property
    def nodes(self) -> NDArray:
        return self.quadrature.measure.nodes

    @property
    def layout(self) -> None:
        return None


@dataclass(frozen=True, eq=False)
class FromBasis:
    r"""The origin of a basis's continuum coefficient axis, on the basis's
    domain, indexed by its layout (ruling 81).

    The kind is the layout's to say, one ``match``: a cell layout indexes
    cell values (``NODAL``, a coordinate cone, unit ones); a degree, a mode,
    a polynomial or a tensor layout indexes spectral coefficients, ``Modal`` with the coefficients of
    the constant function :math:`1` as its unit (ruling 37).
    """

    basis: Basis | TruncatedBasis

    @property
    def kind(self) -> AxisKind:
        match self.basis.layout:
            case CellLayout():
                return BasisKind.NODAL
            case HarmonicLayout() | LegendreLayout() | ModeLayout() | PolynomialLayout() | TensorLayout() as coefficients:
                return Modal(coefficients.unit)
            case layout:
                raise TypeError(
                    f"FromBasis: {type(self.basis).__name__} carries "
                    f"{type(layout).__name__}, which is not a coefficient "
                    f"layout (orpheus.numerics.coefficient_layout)"
                )

    @property
    def support(self) -> Manifold:
        return self.basis.domain

    @property
    def nodes(self) -> None:
        return None

    @property
    def layout(self) -> CoefficientLayout:
        return self.basis.layout


@dataclass(frozen=True, eq=False)
class FromFrame:
    r"""The origin of a moment head a frame has dressed with its Parseval metric.

    The head holds COVARIANT moments :math:`\varphi = \Phi^{\mathsf T} W
    \psi`, so the coordinates of the constant function there are its
    moments, :math:`u = \Phi^{\mathsf T} w`, the trial basis's analysis of
    :math:`1` (ruling 37), not the basis's coefficient :math:`e_0`. The
    frame is the Stage-2 generator: the pairing basis ⊗ measure induces the
    head's metric, and only the frame can re-dress the head at another order
    (:meth:`~orpheus.numerics.frame.GalerkinFrame.at_order`).
    """

    frame: FrameBase

    @property
    def kind(self) -> AxisKind:
        f = self.frame
        match f.basis.layout:
            case ModeLayout() as modes:
                # a mode basis DECLARES its unit (a computed Phi^T w is
                # round-off where the span is orthogonal to the constants)
                return Modal(modes.unit)
            case _:
                return Modal(np.asarray(f.basis.analyze(np.ones(f.measure.n_points), f.table, f.measure.weights)))

    @property
    def support(self) -> Manifold:
        return self.frame.basis.domain

    @property
    def nodes(self) -> None:
        return None

    @property
    def layout(self) -> CoefficientLayout:
        return self.frame.basis.layout


#: What minted an axis, the closed sum of its five arms (ruling 64).
Origin = FromMeasure | FromQuadrature | FromBasis | FromFrame | Declared


@dataclass(frozen=True, eq=False)
class Axis(ContentIdentity):
    r"""One tensor factor of a function space — a frozen value object.

    Parameters
    ----------
    label : str
        The factor's role name (``"energy"``, ``"spatial"``, …). Part of
        the identity: two same-shaped factors with different roles are
        different axes.
    shape : tuple[int, ...]
        The factor's index shape, rank ≥ 1 (refused otherwise).
    weights : NDArray | None, default None
        The factor measure over exactly ``shape``. ``None`` IS the
        counting measure (identity metric) — deliberately, always; an
        all-ones array is CANONICALIZED to ``None`` at construction so
        the counting measure has one spelling and one identity. Stored
        canonicalized (``w + 0.0``, killing ``-0.0``), defensively
        copied, read-only. Non-finite entries are refused; signed
        weights are legal (quadrature families need them).
    origin : Origin
        Keyword-only, no default: what minted the axis (the module
        docstring's table). The kind, the support and the atoms are read
        from it, so the basis character is spelled at every mint and
        cannot disagree with the generator. NOT part of the identity.
        Prefer the generator mints (``Axis.of_measure(measure, label)`` /
        ``Axis.of_quadrature(quad)``); spell ``origin=Declared(kind)`` only
        where no generating object exists.

    Notes
    -----
    Frozen and hashable with **structural equality per subclass**:
    ``__eq__``/``__hash__`` read the content digest of ``(type, label,
    shape, weights, kind, support, nodes)``
    (:class:`~orpheus.values.content.ContentIdentity`); subclasses extend
    the content by declaring fields. The class check keeps an
    ``EnergyAxis`` and a field-identical generic ``Axis`` unequal.
    ``dataclasses.replace(axis, origin=...)`` re-derives the kind, support
    and atoms, so a re-dressing cannot keep a stale kind.
    """

    label: str
    shape: tuple[int, ...]
    weights: NDArray | None = field(default=None, repr=False)
    # Provenance is not identity: ``compare=False`` is the content encoder's
    # one spelling of "not content" (#405 P1 step 5).
    origin: Origin = field(kw_only=True, repr=False, compare=False)
    kind: AxisKind = field(init=False)
    support: Manifold | None = field(init=False, repr=False)
    nodes: NDArray | None = field(init=False, repr=False)
    layout: CoefficientLayout | None = field(init=False, repr=False)

    def __post_init__(self) -> None:
        shape = tuple(int(n) for n in self.shape)
        object.__setattr__(self, "shape", shape)
        if not isinstance(self.origin, Origin):
            raise TypeError(
                f"axis '{self.label}': origin must be one of FromMeasure, "
                f"FromQuadrature, FromBasis, FromFrame, Declared; got {self.origin!r}"
            )
        kind = self.origin.kind
        kind.unit_over(shape)  # a modal unit over another shape is refused
        object.__setattr__(self, "kind", kind)
        object.__setattr__(self, "support", self.origin.support)
        layout = self.origin.layout
        if layout is not None and layout.shape != shape:
            raise ValueError(
                f"axis '{self.label}': the basis's layout indexes shape "
                f"{layout.shape}, not the axis shape {shape}"
            )
        object.__setattr__(self, "layout", layout)
        if len(shape) < 1:
            raise ValueError(
                f"an Axis needs rank >= 1, got shape {shape!r} — a rank-0 "
                f"factor has no index set to measure"
            )
        nodes = self.origin.nodes
        if nodes is not None and np.shape(nodes)[0] != int(np.prod(shape)):
            raise ValueError(
                f"axis '{self.label}': a sampled axis has one atom per index; "
                f"{np.shape(nodes)[0]} atoms against shape {shape}"
            )
        object.__setattr__(self, "nodes", None if nodes is None else _read_only(nodes))
        if self.weights is not None:
            w = np.ascontiguousarray(self.weights, dtype=float)
            if w.shape != shape:
                raise ValueError(
                    f"axis weights must live over exactly the axis shape: "
                    f"weights shape {w.shape} != axis shape {shape}"
                )
            if not bool(np.isfinite(w).all()):
                raise ValueError(
                    f"axis weights must be finite (a factor measure has "
                    f"finite weights); got {w!r}"
                )
            if bool((w == 1.0).all()):
                # Canonicalization: the counting measure has ONE spelling.
                object.__setattr__(self, "weights", None)
            else:
                # ``+ 0.0`` canonicalizes -0.0 -> +0.0 AND forces a fresh
                # allocation (defensive copy even when the input was
                # already a contiguous float array).
                object.__setattr__(self, "weights", _read_only(w))

    @property
    def unit(self) -> NDArray:
        r"""The coordinates of the constant function :math:`1` on this factor,
        read from the kind: ones on a ``NODAL`` axis, the carried unit on a
        :class:`Modal` one."""
        return self.kind.unit_over(self.shape)

    @classmethod
    def of_measure(cls, measure: DiscreteMeasure, label: str) -> Axis:
        r"""The space-factor axis a discrete measure samples.

        A nodal axis is the image of a discrete measure: its weights are the
        measure's, and its content keeps the support and the atoms (ruling
        71), so two measures with equal weights and different atoms, or the
        same atoms on two manifolds, mint two axes. The measure itself rides
        in :class:`FromMeasure` as provenance.

        ``kind`` is NOT a parameter: a discrete measure's samples are point
        values with a coordinate cone, so the origin implies
        :attr:`BasisKind.NODAL`. This reads only the measure's attributes,
        so the module takes no runtime import of
        :mod:`orpheus.numerics.measure`.
        """
        return Axis(label, (measure.n_points,), weights=measure.weights, origin=FromMeasure(measure))

    @classmethod
    def of_quadrature(cls, quadrature: Quadrature, label: str = "angular") -> Axis:
        r"""The angular space-factor axis a directional quadrature samples,
        with THAT rule as its origin.

        The label DEFAULTS: a directional quadrature is the angular
        generator, so its axis's role is intrinsic, and defaulting the
        spelling makes a label-twin across mint sites unspellable (the
        P4-remainder's second-mint hazard). Pass a label only for a
        deliberately re-labelled factor.

        The content is the rule's measure's, so this is the axis
        ``of_measure(quadrature.measure, label)`` mints; the origin is the
        quadrature because the measure alone does not answer everything the
        axis forgot: ``level_indices`` lives on the
        :class:`~orpheus.numerics.quadrature.directional.LevelStructure`
        side-channel, carried by the quadrature. A consumer recovering the
        angular geometry through the space matches
        ``case FromQuadrature(quadrature)`` and gets the whole rule, not a
        projection of it ([M] 2026-08-29, CS5 ground memo section L).
        """
        measure = quadrature.measure
        return Axis(label, (measure.n_points,), weights=measure.weights, origin=FromQuadrature(quadrature))


    @classmethod
    def of_basis(cls, basis: Basis | TruncatedBasis) -> Axis:
        r"""The coefficient axis a basis spans: indexed by its layout, measured
        by its continuum Gram (``None`` for a measure-free basis), on its domain
        (ruling 81). The one generic coefficient axis: the family is the
        layout, the spent axis of a Legendre basis is the domain
        :math:`S^2/O(2)_a`, and no family subclass exists. A frame re-dresses
        it with the discrete Parseval inverse
        (:func:`~orpheus.numerics.frame_faces.basis_space`)."""
        layout = basis.layout
        return Axis(layout.label, layout.shape, weights=basis.coefficient_weights, origin=FromBasis(basis))


def quadrature_refusal(axis: Axis, consumer: str) -> ValueError:
    r"""The refusal a consumer raises when it needs the generating quadrature
    and ``axis`` was not minted through one.

    The consumers that recover the angular geometry through the space (the
    direction cosines, the level fibration) match ``case
    FromQuadrature(quadrature)`` themselves; this is the one spelling of
    their other arm, so the message names both the axis and the asker and a
    refusal is diagnosable from the text alone (P4-remainder G5).
    """
    return ValueError(
        f"axis '{axis.label}': {consumer} needs the generating Quadrature, "
        f"but this axis was minted through {type(axis.origin).__name__}, not "
        f"FromQuadrature. Mint the axis via Axis.of_quadrature(quad) so "
        f"consumers can recover what the axis forgot."
    )


@dataclass(frozen=True, eq=False)
class EnergyAxis(Axis):
    r"""The multigroup energy factor — a 1-D mesh in energy.

    **The faces reading.** A multigroup group structure (the data layer's
    ``EnergyGrid``) is a one-dimensional MESH in energy: the group boundaries are its
    FACES, the groups are its CELLS, and condensation is the mesh-overlap
    map (``EnergyGrid.overlap_to`` — fractional re-binning, fine → coarse
    only). The one-group member is the one-CELL energy mesh — its edges
    and weighting spectrum survive because they define :math:`\bar\sigma`,
    which is exactly what the Bateman/depletion pairing
    :math:`\langle\bar\sigma, \phi\rangle` consumes. The axis therefore
    persists down to its terminal one-cell member (collapse doctrine,
    clause 2: partition-integration of an L¹ field class).

    **The counting-measure theorem.** Multigroup flux components are group
    INTEGRALS (covariant, extensive): :math:`\phi_g = \int_g \phi(E)\,dE`.
    Cross sections are flux-weighted group AVERAGES (contravariant,
    intensive). The convention is chosen so that
    :math:`\int \sigma(E)\,\phi(E)\,dE = \sum_g \sigma_g\,\phi_g`
    EXACTLY — no group widths appear. The energy metric is therefore the
    COUNTING measure **as a theorem, not a default**: metric = I,
    :math:`V \cong V^*` isometrically, and the adjoint along energy is
    the plain transpose. Construction enforces it: a weighted
    ``EnergyAxis`` is refused (use a generic :class:`Axis` for
    deliberately non-physical toys), and both constructors mint
    ``weights=None``.

    **The V/V\* collapse hook (declared now, built at S7/Campaign 2).**
    Condensation acts as plain SUM on V (integrals add) and as
    flux-weighted AVERAGE on V* (averages re-weight); collapse
    adjoint-consistency is precisely that pair being mutually adjoint
    under the counting pairing. This axis records the group structure
    those morphisms will consume.

    **Identity (ruled Q2): ng + edges CONTENT.** ``of_edges`` axes carry
    the boundary energies; equality reads their BYTES, never the identity
    of whatever object produced them (a group structure in the data layer
    is an ``eq=False`` value minted fresh per access, so two mints from
    one producer must, and do, yield equal axes). ``synthetic(ng)`` axes
    carry no edges; identity is ``ng`` alone, and a synthetic axis NEVER
    equals an ``of_edges`` axis of the same ``ng`` (same index set,
    different partition — different axis).

    Edges follow the canonical fast-first convention (strictly
    DESCENDING, group 0 = fastest). The invariant is the data layer's
    group-structure contract, checked there once, at construction, and
    not re-checked here; the production route to a gridded axis is the
    grid's own view, ``EnergyGrid.as_axis``.
    """

    edges: NDArray | None = field(default=None, kw_only=True, repr=False)
    #: Fixed, not a parameter: groups are the CELLS of a 1-D energy mesh, so
    #: the factor is NODAL with no atoms (its edges are its content; ruling 82).
    origin: Origin = field(default=Declared(BasisKind.NODAL), init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        super().__post_init__()
        if self.weights is not None:
            raise ValueError(
                "an EnergyAxis cannot carry weights: the multigroup "
                "convention makes the energy measure COUNTING as a theorem "
                "(group integrals x group averages pair without widths). "
                "For a deliberately non-physical weighted toy, use a "
                "generic Axis."
            )
        if len(self.shape) != 1:
            raise ValueError(
                f"an EnergyAxis is rank 1 (groups index a 1-D energy "
                f"mesh); got shape {self.shape}"
            )
        if self.edges is not None:
            edges = np.ascontiguousarray(self.edges, dtype=float) + 0.0
            edges.setflags(write=False)
            object.__setattr__(self, "edges", edges)
            (ng,) = self.shape
            if edges.shape != (ng + 1,):
                raise ValueError(
                    f"EnergyAxis edges must be the {ng + 1} group "
                    f"boundaries of its {ng} groups; got shape "
                    f"{edges.shape}"
                )

    @classmethod
    def of_edges(cls, edges: NDArray) -> "EnergyAxis":
        """The axis of a real group structure, minted from its boundary
        energies — identity = ng + edges bytes.

        The ``ng + 1`` edges are the axis's CONTENT, not a generating
        object, so the origin stays the class's fixed ``Declared`` arm; the
        group count is derived (``len(edges) - 1``) and the label is fixed
        at ``"energy"``, so neither a count mismatch nor a label twin is
        spellable. Monotonicity is NOT re-checked: it is the group
        structure's contract, checked once where that structure is built
        (the production route is the grid's view, ``EnergyGrid.as_axis``).
        """
        edges = np.asarray(edges, dtype=float)
        if edges.ndim != 1 or edges.shape[0] < 2:
            raise ValueError(
                f"EnergyAxis.of_edges needs a 1-D array of >= 2 group "
                f"boundaries; got shape {edges.shape}"
            )
        return cls("energy", (edges.shape[0] - 1,), edges=edges)

    @classmethod
    def synthetic(cls, ng: int) -> "EnergyAxis":
        """The axis of a grid-less ``ng``-group problem — identity = ng only.

        The honest spelling for fixtures/libraries that declare a group
        COUNT with no boundary energies (every shipped ``get_mixture``
        pair has ``eg is None``). Deliberately UNEQUAL to any
        ``of_edges`` axis of the same ``ng``: same index set, no
        partition data — a different axis.
        """
        return cls("energy", (int(ng),), edges=None)
