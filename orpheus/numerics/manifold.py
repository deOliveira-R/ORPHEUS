r"""The manifold a measure lives on and a basis is defined over.

A :class:`Manifold` is the *point set* — level 1 of the three-level stack
that :doc:`the theory page </theory/foundations/spaces>` describes:

.. code-block:: text

    level 1   the manifold M            S^2, [-1,1], R^d, M/H, M x N
    level 2   fields on it              L^2(M)          <- FunctionSpace
    level 3   coefficients              R^K             <- a Basis's shape

Until 2026-08 level 1 was ``Space = str`` — an opaque tag whose own
docstring called the entries *"recommendations, not constraints"*.  Three
consequences, all measured, and all of them the reason this type exists:

* **Nothing could be refused.**  A 1-D quadrature could declare its nodes
  live on :math:`S^2` while supplying :math:`\mu_y = \mu_z = 0`, which is
  not a point of the sphere.  That forgery is the first link of
  :ref:`ERR-080 <vv-error-ERR-080>` — a live wrong-answer defect in
  :math:`P_L` scattering.  :meth:`Manifold.contains` refuses it at
  construction, three hops before the symptom.
* **The algebra ran as string concatenation.**  ``f"{a} x {b}"`` was the
  product, ``f"{a}/{H.name}"`` the quotient, ``f"phi_*({a})"`` the
  pushforward codomain.  Those are :meth:`__mul__`, :meth:`quotient` and a
  chart's codomain; the interpolation *was* the operation, unnamed.
* **The vocabulary drifted.**  ``'S^2/<sigma_y>'`` and ``'S^2/sigma_y'``
  both shipped, naming one quotient and unequal under ``==``.

The type is a **closed sum**: a small frozen variant set, with the
operations every manifold answers (:attr:`dim`, :meth:`contains`,
:meth:`__mul__`, :attr:`name`) on the base, and the operations only a
*quotient* can answer on :class:`Quotient` alone.  A sphere has no syzygy
ideal, and under this shape it cannot be asked for one.

⚠ **This module imports nothing from** :mod:`orpheus.numerics` **at
runtime, and that is load-bearing, not tidiness.**
:mod:`orpheus.numerics.symmetry` imports :mod:`orpheus.numerics.measure`
at module scope (``symmetry.py:98``), and ``measure`` will import this
module — so a module-scope ``manifold -> symmetry`` edge closes the cycle
``measure -> manifold -> symmetry -> measure`` and breaks
``import orpheus.numerics.measure``.  :class:`SubgroupOfO3` is therefore
referenced under :data:`typing.TYPE_CHECKING` only; nothing here calls a
method on it.  (``measure.py:1005`` defers its own ``symmetry`` import to
function scope for exactly the same reason.)

Why a closed sum rather than a polymorphic hierarchy: the members are
stable (about a dozen, and every new orbit space is another *instance* of
:class:`Quotient`, not another class) while the operations keep growing.
The deciding evidence is that all of the derivation fields a catalogued
quotient carries belong to :class:`Quotient` alone — on a hierarchy they
would sit on the base returning ``None`` for every non-quotient, which is
the tax :attr:`SubgroupOfO3.mirror_axis` already pays and which
``directional.py`` already branches on.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

if TYPE_CHECKING:  # runtime import would close the cycle documented above
    from orpheus.numerics.symmetry import SubgroupOfO3

__all__ = [
    "Manifold",
    "Sphere",
    "Circle",
    "Interval",
    "RealSpace",
    "EnergyGroups",
    "IndexSet",
    "Product",
    "Quotient",
    "SPHERE",
    "CIRCLE",
    "UNIT_INTERVAL",
    "COSINE_INTERVAL",
    "HALF_LINE",
    "REAL_LINE",
    "ENERGY",
]

#: Tolerance for :meth:`Manifold.contains` on a curved manifold.  A node
#: set is admitted when it satisfies the defining equation to this
#: absolute tolerance.  Chosen to match the construction quality of the
#: shipped quadrature rules, whose nodes are exact to a few ULP; it is
#: NOT a physics tolerance and no caller should widen it to make a
#: measure fit.
_MEMBERSHIP_ATOL = 1e-12


@dataclass(frozen=True)
class Manifold(ABC):
    r"""The point set a measure is supported on — the closed sum's base.

    Only the operations **every** manifold answers live here.  Anything a
    particular family alone can answer belongs on that family, so that
    asking a sphere for a quotient's syzygy ideal is a type error rather
    than a ``None``.

    Subclassing outside this module is not supported: the set is closed
    so that a new operation can be checked against every member.
    """

    @property
    @abstractmethod
    def dim(self) -> int:
        """Topological dimension of the manifold."""
        ...

    @property
    @abstractmethod
    def name(self) -> str:
        """Display name, used in messages and in ``L2[...]`` space names.

        The names reproduce the retired ``SPACE_*`` string tags exactly,
        so the migration off ``Space = str`` does not silently re-word
        any message a test pins.
        """
        ...

    @abstractmethod
    def contains(self, points: ArrayLike) -> bool:
        r"""``True`` iff every row of ``points`` is a point OF this manifold.

        This is the predicate that makes a support claim falsifiable.  A
        measure whose nodes fail it is not a measure on this manifold,
        whatever its tag says.

        Parameters
        ----------
        points
            ``(n, d)`` array of candidate points, or ``(n,)`` for a
            1-dimensional ambient space.

        Notes
        -----
        The check is on the AMBIENT coordinates, so it catches the
        ERR-080 class directly: a 1-D quadrature's ordinates
        :math:`(\mu, 0, 0)` satisfy :math:`\|\Omega\| = 1` only at
        :math:`\mu = \pm 1`, so :class:`Sphere` refuses them.
        """
        ...

    def __mul__(self, other: object) -> "Product":
        r"""The product manifold :math:`M \times N`.

        Replaces ``new_space = f"{self.support} x {other.support}"``.
        """
        if not isinstance(other, Manifold):
            return NotImplemented
        return Product(self, other)

    def quotient(self, group: "SubgroupOfO3") -> "Quotient":
        r"""The orbit space :math:`M/H`, by catalogue lookup.

        Computing an orbit space from scratch — invariant ring, syzygy
        ideal by elimination, then the Procesi–Schwarz PSD condition — is
        mechanical but is a symbolic-algebra engine, and the groups that
        occur in transport number about a dozen.  So entries are derived
        once by hand and recorded; the engine is **deferred, not
        refused**, and when it is built it populates exactly these
        fields rather than introducing a parallel representation.

        Raises
        ------
        NotImplementedError
            When no entry exists, naming the missing WORK rather than
            the gap, so a fresh session can pick it up.
        """
        entry = _ORBIT_CATALOGUE.get((type(self), group.name))
        if entry is None:
            raise NotImplementedError(
                f"no catalogue entry for {self.name}/{group.name}: derive "
                f"it once (minimal invariants p_1..p_k of R[x]^H; syzygy "
                f"ideal I by elimination; Procesi-Schwarz P_ij = "
                f"<grad p_i, grad p_j> with P >= 0; intersect with the "
                f"ideal of {self.name}) and register it in "
                f"orpheus/numerics/manifold.py's _ORBIT_CATALOGUE, or "
                f"implement the derivation engine. "
                f"Catalogued today (manifold CLASS / group): "
                f"{sorted(f'{m.__name__}/{g}' for m, g in _ORBIT_CATALOGUE)}."
            )
        return entry(self, group)

    def _as_points(self, points: ArrayLike, ambient: int) -> NDArray:
        """Coerce ``points`` to ``(n, ambient)``, or raise saying why not."""
        arr = np.atleast_1d(np.asarray(points, dtype=float))
        if ambient == 1 and arr.ndim == 1:
            arr = arr[:, None]
        if arr.ndim != 2 or arr.shape[1] != ambient:
            raise ValueError(
                f"{self.name}.contains expects points of ambient dimension "
                f"{ambient} — an (n, {ambient}) array"
                f"{' or an (n,) array' if ambient == 1 else ''}; got shape "
                f"{arr.shape}."
            )
        return arr


# ---------------------------------------------------------------------------
# Curved members — the defining equation IS the membership predicate
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Sphere(Manifold):
    r""":math:`S^2 = \{\Omega \in \mathbb{R}^3 : \|\Omega\| = 1\}`."""

    @property
    def dim(self) -> int:
        return 2

    @property
    def name(self) -> str:
        return "S^2"

    def contains(self, points: ArrayLike) -> bool:
        arr = self._as_points(points, 3)
        return bool(
            np.all(
                np.abs(np.linalg.norm(arr, axis=1) - 1.0) <= _MEMBERSHIP_ATOL
            )
        )


@dataclass(frozen=True)
class Circle(Manifold):
    r""":math:`S^1 = \{u \in \mathbb{R}^2 : \|u\| = 1\}`.

    The circle names the MANIFOLD, not a chart of it — the tag was
    ``"[0,2pi)"`` until 2026-08-02, which asserted a coordinate.  Its
    nodes are the roots of unity as POINTS, because only in that
    representation is an on-axis node exactly ``0.0``.
    """

    @property
    def dim(self) -> int:
        return 1

    @property
    def name(self) -> str:
        return "S^1"

    def contains(self, points: ArrayLike) -> bool:
        arr = self._as_points(points, 2)
        return bool(
            np.all(
                np.abs(np.linalg.norm(arr, axis=1) - 1.0) <= _MEMBERSHIP_ATOL
            )
        )


# ---------------------------------------------------------------------------
# Flat members
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Interval(Manifold):
    r"""The real interval :math:`[a, b]`, endpoints included.

    One family, not three tags: the retired ``SPACE_INTERVAL_M11``,
    ``SPACE_INTERVAL_01``, ``SPACE_HALF_LINE`` and ``SPACE_R`` are four
    MEMBERS of it — which is what ``generating_measure.py``'s
    ``support=f"[{a},{b}]"`` was already saying, one interpolation at a
    time.
    """

    a: float
    b: float

    def __post_init__(self) -> None:
        if not self.a < self.b:
            raise ValueError(
                f"Interval requires a < b; got a={self.a!r}, b={self.b!r}. "
                f"A degenerate interval is a point, not a 1-manifold."
            )

    @property
    def dim(self) -> int:
        return 1

    @property
    def name(self) -> str:
        # Reproduce the retired tags verbatim so no pinned message moves.
        if self.a == -np.inf and self.b == np.inf:
            return "R"
        if self.a == 0.0 and self.b == np.inf:
            return "[0,inf)"
        return f"[{_fmt(self.a)},{_fmt(self.b)}]"

    def contains(self, points: ArrayLike) -> bool:
        arr = self._as_points(points, 1)[:, 0]
        # Finiteness is a SEPARATE condition from the bounds: on an
        # unbounded interval `inf <= inf` holds, and infinity is not a
        # point of R.  A NaN fails every comparison and so is refused
        # by the bounds anyway; asserting finiteness says why.
        return bool(
            np.all(np.isfinite(arr))
            and np.all(arr >= self.a - _MEMBERSHIP_ATOL)
            and np.all(arr <= self.b + _MEMBERSHIP_ATOL)
        )


@dataclass(frozen=True)
class RealSpace(Manifold):
    r""":math:`\mathbb{R}^d` — the spatial manifold, indexed by ``d``.

    Replaces the interpolated ``support=f"spatial_R{ndim}"``.
    """

    d: int

    def __post_init__(self) -> None:
        if self.d < 1:
            raise ValueError(
                f"RealSpace requires d >= 1; got d={self.d!r}."
            )

    @property
    def dim(self) -> int:
        return self.d

    @property
    def name(self) -> str:
        return f"spatial_R{self.d}"

    def contains(self, points: ArrayLike) -> bool:
        arr = self._as_points(points, self.d)
        return bool(np.all(np.isfinite(arr)))


# ---------------------------------------------------------------------------
# Discrete members — a finite index set carries no geometry
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class IndexSet(Manifold):
    r"""A finite index set — ``n`` points with no metric structure.

    Two shipped sites minted this concept under incompatible spellings:
    ``frame.py``'s ``f"index({axis_label})"`` and
    ``loss_kernel_gauge.py``'s ``f"sn_trace_orbit{orbit}_g{group}"``,
    whose "points" are trace DOF indices cast to float.  They are one
    family; ``label`` is what distinguishes the instances.

    ``n`` is optional: a label alone identifies the set, and the
    cardinality is often the measure's own ``n_points``.
    """

    label: str
    n: int | None = None

    @property
    def dim(self) -> int:
        return 0

    @property
    def name(self) -> str:
        return f"index({self.label})"

    def contains(self, points: ArrayLike) -> bool:
        arr = self._as_points(points, 1)[:, 0]
        if not np.all(arr == np.floor(arr)):
            return False
        if self.n is None:
            return bool(np.all(arr >= 0))
        return bool(np.all((arr >= 0) & (arr < self.n)))


@dataclass(frozen=True)
class EnergyGroups(Manifold):
    r"""The multigroup energy axis — a counting set of ``ng`` groups.

    Distinct from a bare :class:`IndexSet` because the measure's PHASE
    depends on it: an energy group index and any other integer-noded
    counting rule are indistinguishable from their nodes alone, which is
    exactly why the tag had to supply the physical identity.
    """

    ng: int | None = None

    @property
    def dim(self) -> int:
        return 0

    @property
    def name(self) -> str:
        return "energy"

    def contains(self, points: ArrayLike) -> bool:
        arr = self._as_points(points, 1)[:, 0]
        if not np.all(arr == np.floor(arr)):
            return False
        if self.ng is None:
            return bool(np.all(arr >= 0))
        return bool(np.all((arr >= 0) & (arr < self.ng)))


# ---------------------------------------------------------------------------
# Recursive members
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Product(Manifold):
    r""":math:`M \times N`, with coordinates concatenated in that order."""

    left: Manifold
    right: Manifold

    @property
    def dim(self) -> int:
        return self.left.dim + self.right.dim

    @property
    def name(self) -> str:
        return f"{self.left.name} × {self.right.name}"

    def contains(self, points: ArrayLike) -> bool:
        arr = np.atleast_2d(np.asarray(points, dtype=float))
        split = _ambient(self.left)
        if arr.shape[1] != split + _ambient(self.right):
            raise ValueError(
                f"{self.name}.contains expects ambient dimension "
                f"{split + _ambient(self.right)}; got shape {arr.shape}."
            )
        return self.left.contains(arr[:, :split]) and self.right.contains(
            arr[:, split:]
        )


@dataclass(frozen=True)
class Quotient(Manifold):
    r"""The orbit space :math:`M/H`, and the derivation that produced it.

    ⭐ **The fields below ARE the derivation procedure's output**, not a
    human summary of its answer.  That is what makes this class the SEED
    of the deferred engine rather than a rival to it: an engine ships
    when it computes these same fields, introducing no new type and no
    translation layer.  The ~12 catalogued entries' symbolic regression
    tests are therefore the engine's specification, written before it.

    ⚠ An orbit space is generally an **orbifold**, not a quotient
    manifold: where the action is not free, the image of :math:`H`'s
    fixed-point set is a singular stratum.  Any design that assumes a
    quotient is a smooth submersion is wrong exactly there — and for
    :math:`S^2/SO(2)` that locus is :math:`\mu = \pm 1`, the poles.
    """

    base: Manifold
    by: "SubgroupOfO3"
    #: The orbit space realized as an honest manifold — what a chart of
    #: ``M/H`` maps ONTO.  For ``S^2/SO(2)`` this is ``Interval(-1, 1)``.
    realization: Manifold
    #: Minimal generators of the invariant ring, as SymPy expressions in
    #: the ambient coordinates.  ``()`` is a legal value only for the
    #: trivial group.
    generators: tuple[Any, ...] = ()
    #: Generators of the syzygy ideal ``I = ker(R[y] -> R[x])``.  Empty
    #: when the invariants are algebraically independent.
    syzygy: tuple[Any, ...] = ()
    #: ``P_ij = <grad p_i, grad p_j>`` re-expressed in the invariants —
    #: the Procesi–Schwarz matrix whose ``P >= 0`` locus cuts out the
    #: orbit space's inequalities.
    gram: Any | None = None
    #: ``det P``.  Its vanishing locus IS the orbit-space boundary, and
    #: it is the squared orbit radius.
    det_gram: Any | None = None
    #: ``"hand"`` for an entry derived once by a human, ``"engine"`` once
    #: the derivation engine populates it.  A MIXED state is expressible
    #: on purpose: an incremental engine rollout is exactly that, and
    #: without this field the migration would have to be all-or-nothing.
    derived_by: str = "hand"

    @property
    def dim(self) -> int:
        return self.realization.dim

    @property
    def name(self) -> str:
        return f"{self.base.name}/{self.by.name}"

    def contains(self, points: ArrayLike) -> bool:
        """Membership is decided in the REALIZATION's coordinates."""
        return self.realization.contains(points)

    @property
    def is_free(self) -> bool:
        """``True`` iff the action has no fixed points, i.e. no stratum."""
        return self.singular_stratum == ()

    #: Points of the realization where ``det P`` vanishes — the singular
    #: stratum, DERIVED rather than declared (nothing declares what it
    #: can derive).  Recorded as the explicit locus for a catalogued
    #: entry; an engine solves ``det P = 0`` to populate it.
    singular_stratum: tuple[float, ...] = ()


# ---------------------------------------------------------------------------
# The shipped members, under their retired tag names
# ---------------------------------------------------------------------------

SPHERE = Sphere()
CIRCLE = Circle()
COSINE_INTERVAL = Interval(-1.0, 1.0)  # was SPACE_INTERVAL_M11
UNIT_INTERVAL = Interval(0.0, 1.0)  # was SPACE_INTERVAL_01
HALF_LINE = Interval(0.0, float(np.inf))  # was SPACE_HALF_LINE
REAL_LINE = Interval(float(-np.inf), float(np.inf))  # was SPACE_R
ENERGY = EnergyGroups()


# ---------------------------------------------------------------------------
# The orbit-space catalogue
# ---------------------------------------------------------------------------


def _sphere_mod_so2(base: Manifold, group: "SubgroupOfO3") -> Quotient:
    r"""``S^2 / SO(2) = [-1, 1]``, derived per the standard procedure.

    Invariants of the axial rotations acting on :math:`\mathbb{R}^3` are
    :math:`p_1 = z` and :math:`p_2 = x^2 + y^2`; they are algebraically
    independent, so the syzygy ideal is empty.  The Procesi–Schwarz
    matrix is

    .. math:: P = \begin{pmatrix} 1 & 0 \\ 0 & 4 p_2 \end{pmatrix},
              \qquad \det P = 4 p_2,

    so :math:`\mathbb{R}^3/SO(2) = \{p_2 \ge 0\}`.  Adjoining the sphere's
    own ideal :math:`p_1^2 + p_2 = 1` and writing :math:`\mu = p_1`:

    .. math:: S^2/SO(2) = \{\mu \in \mathbb{R} : 1 - \mu^2 \ge 0\}
              = [-1, 1].

    ⭐ :math:`\det P = 4(1-\mu^2)` is the squared orbit radius, and it is
    the SAME polynomial as the curvilinear angular-redistribution
    coefficient in :math:`(1/r)\,\partial_\mu[(1-\mu^2)\psi]` — the
    quotient's connection term.  Its zeros, :math:`\mu = \pm 1`, are the
    poles: the two points with full stabilizer, hence the singular
    stratum, hence why the orbit space is an orbifold.

    SymPy is imported lazily: it costs ~250 ms and nothing else in
    :mod:`orpheus.numerics` needs it, so a session that never asks for a
    quotient never pays for it.
    """
    import sympy as sp

    x, y, z = sp.symbols("x y z", real=True)
    p1, p2 = sp.symbols("p1 p2", real=True)

    invariants = (z, x**2 + y**2)
    grad = [[sp.diff(p, v) for v in (x, y, z)] for p in invariants]
    gram_xyz = sp.Matrix(
        [[sum(gi * gj for gi, gj in zip(a, b)) for b in grad] for a in grad]
    )
    # Re-express in the invariants: <grad p2, grad p2> = 4(x^2+y^2) = 4 p2.
    gram = gram_xyz.subs({x**2 + y**2: p2}).applyfunc(sp.simplify)
    gram = gram.subs({x: sp.sqrt(p2), y: 0})  # on the orbit, wlog
    gram = gram.applyfunc(sp.simplify)

    return Quotient(
        base=base,
        by=group,
        realization=COSINE_INTERVAL,
        generators=(p1 - z, p2 - (x**2 + y**2)),
        syzygy=(),
        gram=gram,
        det_gram=sp.simplify(gram.det()),
        derived_by="hand",
        singular_stratum=(-1.0, 1.0),
    )


#: ``(manifold type, group name) -> builder``.  Keyed on the PAIR because
#: a quotient is binary dispatch: it is a property of neither operand
#: alone.  About a dozen entries are expected; the engine that would
#: compute them instead of reading them is deferred, not refused.
_ORBIT_CATALOGUE: dict[tuple[type, str], Any] = {
    (Sphere, "SO2"): _sphere_mod_so2,
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _ambient(m: Manifold) -> int:
    """Ambient coordinate count — how many columns ``m`` consumes."""
    match m:
        case Sphere():
            return 3
        case Circle():
            return 2
        case Interval() | IndexSet() | EnergyGroups():
            return 1
        case RealSpace(d=d):
            return d
        case Product(left=left, right=right):
            return _ambient(left) + _ambient(right)
        case Quotient(realization=realization):
            return _ambient(realization)
    raise NotImplementedError(
        f"ambient dimension is undefined for {type(m).__name__}. A new "
        f"Manifold variant must be added here — this match is deliberately "
        f"exhaustive so that adding a member cannot silently skip it."
    )


def _fmt(v: float) -> str:
    """Format an endpoint the way the retired string tags did."""
    return str(int(v)) if float(v).is_integer() else repr(v)
