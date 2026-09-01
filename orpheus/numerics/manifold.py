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
  **ERR-080** — a live wrong-answer defect in :math:`P_L` scattering.
  (Plain text, not a ``:ref:``: `[M]` the ``error-entry`` directive emits a
  ``container`` + ``rubric`` and **no** ``nodes.target``, so there is no
  label to resolve — ``vv-error-ERR-080`` occurs nowhere else in ``docs/``
  or ``orpheus/``.  Every one of the ~150 ERR citations in ``orpheus/`` is
  plain text for the same reason.)  :meth:`Manifold.contains` refuses it at
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
        # M/{e} = M is a THEOREM for every manifold, not a table row —
        # so it is derived here rather than needing one entry per member.
        # Run the same procedure on the trivial group: its invariant ring
        # is the whole polynomial ring, generated by the coordinates
        # themselves, so P_ij = <grad x_i, grad x_j> = delta_ij, and
        # det P = 1 vanishes NOWHERE.  No vanishing locus means no
        # singular stratum, which means the action is free — and the
        # trivial action is free, vacuously, since its only element is
        # the identity.  The procedure reproduces the known answer, so
        # this case is also a positive control on the machinery.
        if group.name == "Trivial":
            return _mod_trivial(self, group)

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
class Ball(Manifold):
    r"""The closed unit ball :math:`D^d = \{p \in \mathbb{R}^d : \|p\| \le 1\}`.

    Minted 2026-08-31 because the derivation demanded it: ``S^2/sigma_y``
    IS the closed 2-disk in invariant coordinates, and no shipped member
    could say so.  ``Product(COSINE_INTERVAL, COSINE_INTERVAL)`` is the
    bounding SQUARE, and the discriminator is measured — ``(0.9, 0.9)``
    is in the square and not in the disk, and it corresponds to NO
    direction, since ``eta^2 + mu^2 = 1.62 > 1`` forces ``xi^2 < 0``.
    """

    d: int

    @property
    def dim(self) -> int:
        return self.d

    @property
    def name(self) -> str:
        return f"D^{self.d}"

    def contains(self, points: ArrayLike) -> bool:
        arr = self._as_points(points, self.d)
        return bool(
            np.all(np.isfinite(arr))
            and np.all(
                np.sum(arr * arr, axis=1) <= 1.0 + _MEMBERSHIP_ATOL
            )
        )


@dataclass(frozen=True)
class FundamentalDomain(Manifold):
    r"""A strict fundamental domain of :math:`H` acting on ``base`` — the
    IMAGE of a section of :math:`M \to M/H`, in the BASE's coordinates.

    Cut from ``base`` by closed half-spaces
    :math:`\langle p, n_i \rangle \ge 0`.  One normal cuts a half-space;
    an antipodal PAIR :math:`\{n, -n\}` cuts the hyperplane
    :math:`\langle p, n \rangle = 0`, so equalities need no second slot —
    which is why one tuple expresses both the :math:`\sigma_y` hemisphere
    and the :math:`SO(2)` half-meridian.

    ⭐ **This is the OTHER half of a quotient, and the tree has always
    needed it.**  :meth:`~orpheus.numerics.measure.DiscreteMeasure.quotient`
    computes orbit representatives and then keeps ``nodes[representative]``
    — a SELECTION, applying no chart — so every measure it emits carries
    the base's ambient columns.  A chart codomain
    (:attr:`Quotient.realization`) cannot validate those points; only a
    fundamental domain can.

    ⛔ **The inequalities must be CLOSED.**  ``[M]`` the cylindrical march
    seeds a level at :math:`\xi = 0` exactly — on the stratum — so a strict
    :math:`\langle p, n \rangle > 0` would refuse a direction production
    marches from (`coding-elegance` anti-pattern #18's half (ii): every
    legal value must be admitted, which is a claim about the PRODUCERS).
    """

    base: Manifold
    #: Outward normals of the closed half-spaces, in the base's ambient
    #: coordinates.  An antipodal pair spells an equality.
    normals: tuple[tuple[float, ...], ...]
    #: Short label for :attr:`name` (e.g. ``"xi>=0"``), since a normal
    #: tuple makes an unreadable one.
    label: str

    @property
    def dim(self) -> int:
        """``base.dim`` less one per antipodal pair (each is an equality).

        An inequality carves a region WITH BOUNDARY and does not drop
        dimension; an equality does.  This is what lets the σ_y
        hemisphere read 2 and the SO(2) half-meridian read 1 from one
        rule — and :class:`Quotient` gates it against the chart.
        """
        seen = {tuple(n) for n in self.normals}
        pairs = sum(
            1
            for n in seen
            if tuple(-c for c in n) in seen and n < tuple(-c for c in n)
        )
        return self.base.dim - pairs

    @property
    def name(self) -> str:
        return f"{self.base.name}|{self.label}"

    def contains(self, points: ArrayLike) -> bool:
        arr = self._as_points(points, _ambient(self.base))
        if not self.base.contains(arr):
            return False
        # Closed half-spaces — see the class docstring's ⛔ note on why
        # this may not be strict.
        return bool(
            all(
                np.all(arr @ np.asarray(n, dtype=float) >= -_MEMBERSHIP_ATOL)
                for n in self.normals
            )
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
    #: Its coordinates are the INVARIANTS', the same language as
    #: :attr:`generators`, :attr:`gram` and :attr:`det_gram`.
    realization: Manifold
    #: The section's IMAGE, in the BASE's coordinates — the other of the
    #: quotient's two honest coordinate systems (user ruling,
    #: 2026-08-31).  ``None`` when no canonical section exists, which for
    #: a positive-dimensional group is the normal case: any half-meridian
    #: sections ``S^2 -> S^2/SO(2)`` and none is distinguished.
    #:
    #: ⭐ Why a SECOND slot rather than a wider ``realization``: the two
    #: answer different questions and the tree needs both.  ``[M]`` the
    #: chart is Mode-12 BLIND to the ERR-080 forgery — ``(x,y,z) -> (x,z)``
    #: drops exactly the coordinate the forgery corrupts, so the forged
    #: row ``(mu, 0)`` is a legal point of the disk — while the section
    #: refuses it.  And ``[M]`` ERR-080 IS a botched section of
    #: ``S^2/SO(2)``: a consumer needed one, the realization is a chart,
    #: and the tree fabricated one by zero-padding to ``(mu, 0, 0)``,
    #: which is off ``S^2``.  With this slot that padding has nowhere to
    #: live.
    fundamental_domain: Manifold | None = None
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

    def __post_init__(self) -> None:
        r"""The two coordinate systems must describe the SAME orbit space.

        A chart codomain and a fundamental domain are two honest views of
        one object, so their dimensions must agree — and that is a real
        check, not a tautology: the fundamental domain derives its ``dim``
        from the base less one per antipodal normal pair, while the
        realization states its own.  ``[M]`` the sigma_y hemisphere reads
        2 against the disk's 2, and an ``SO(2)`` half-meridian would read
        1 against ``[-1,1]``'s 1 — so a mis-specified entry (a hemisphere
        offered for a 1-D orbit space, an equality normal forgotten) is
        refused where it is written, not where it is read.
        """
        if self.fundamental_domain is None:
            return
        if self.fundamental_domain.dim != self.realization.dim:
            raise ValueError(
                f"{self.name}: the fundamental domain "
                f"{self.fundamental_domain.name!r} has dim "
                f"{self.fundamental_domain.dim} but the realization "
                f"{self.realization.name!r} has dim "
                f"{self.realization.dim} — the two must describe the same "
                f"orbit space. Check the normals: an antipodal PAIR spells "
                f"an equality and drops a dimension; a lone normal does not."
            )

    @property
    def dim(self) -> int:
        return self.realization.dim

    @property
    def name(self) -> str:
        return f"{self.base.name}/{self.by.name}"

    def contains(self, points: ArrayLike) -> bool:
        r"""Membership, in EITHER of the quotient's two coordinate systems.

        A point of :math:`M/H` may be given as chart coordinates (the
        :attr:`realization`'s language) or as a representative in the base
        (the :attr:`fundamental_domain`'s).  Both are honest and the type
        knows both, so it accepts both and dispatches on the ambient
        width — the one place the distinction is a genuine local split
        rather than a repeated tag test.

        ⚠ :func:`_ambient` still reports the REALIZATION's width: that is
        the canonical coordinate for composition (a :class:`Product`
        factor must have one width).  This method is deliberately the
        wider of the two.
        """
        arr = np.atleast_1d(np.asarray(points, dtype=float))
        chart = _ambient(self.realization)
        if self.fundamental_domain is not None:
            section = _ambient(self.fundamental_domain)
            width = arr.shape[1] if arr.ndim == 2 else 1
            if width == section and section != chart:
                return self.fundamental_domain.contains(arr)
        return self.realization.contains(arr)

    @property
    def is_free(self) -> bool:
        """``True`` iff the action has no fixed points, i.e. no stratum."""
        return self.singular_stratum is None

    #: The singular stratum as a LOCUS: a SymPy expression in the
    #: realization's coordinates whose vanishing set is the stratum, or
    #: ``None`` for a free action.
    #:
    #: ⛔ This was ``tuple[float, ...]`` until 2026-08-31 and could not
    #: hold the second catalogued entry: ``S^2/sigma_y``'s stratum is the
    #: disk's boundary CIRCLE, not a finite point set.  The first entry's
    #: shape had become the field's type — a stratum is a locus, and two
    #: poles are a locus that happens to be finite.
    #:
    #: It is derivation OUTPUT, not a stored copy of :attr:`det_gram`:
    #: recovering it needs the BASE's own ideal (``det P = 4 p_2`` becomes
    #: ``4(1 - mu^2)`` only after substituting ``p_1^2 + p_2 = 1``), and a
    #: :class:`Quotient` does not carry that ideal.  So the type cannot
    #: recompute it, which is exactly when storing is right.
    singular_stratum: Any | None = None


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
        # No CANONICAL section: any half-meridian sections S^2 -> S^2/SO(2)
        # and none is distinguished, which is the normal case for a
        # positive-dimensional group.  ⛔ ERR-080 is the tree having needed
        # one anyway and fabricated it by zero-padding to (mu, 0, 0), which
        # is off S^2.  Declaring one (the phi=0 half-meridian
        # mu -> (sqrt(1-mu^2), 0, mu)) is a CHOICE and belongs to the step
        # that makes the slab declare its quotient, not to this derivation.
        fundamental_domain=None,
        # det P = 4 p_2 = 4(1 - mu^2) after the sphere ideal p_1^2 + p_2 = 1;
        # it vanishes at the poles mu = +-1.  A locus, in the realization's
        # coordinate — the two poles are a stratum that happens to be finite.
        singular_stratum=sp.simplify(1 - sp.Symbol("u0", real=True) ** 2),
    )


def _sphere_mod_mirror(base: Manifold, group: "SubgroupOfO3") -> Quotient:
    r"""``S^2 / <sigma_a> = D^2``, derived per the standard procedure.

    The shipped CYLINDRICAL FOLD.  ``[M]``
    :meth:`~orpheus.numerics.quadrature.directional.Quadrature.folded_product`
    ends in ``full.quotient(SubgroupOfO3.Mirror("y"))``, so this entry is
    what makes the fold's support expressible as a typed quotient at all.

    Write :math:`a` for the mirrored axis and :math:`b, c` for the other
    two.  The invariant ring :math:`\mathbb{R}[x,y,z]^{\langle \sigma_a
    \rangle}` is minimally generated by :math:`p_1 = x_b`, :math:`p_2 =
    x_c`, :math:`p_3 = x_a^2`.  ⭐ The syzygy ideal is **empty**, and
    predictably so rather than by luck: :math:`\sigma_a` is a REFLECTION,
    so Chevalley–Shephard–Todd forces the invariant ring to be a
    polynomial ring — the three generators are algebraically independent
    by a theorem, not by inspection.  Hence

    .. math:: P = \mathrm{diag}(1,\, 1,\, 4 p_3),
              \qquad \det P = 4 p_3,

    so :math:`\mathbb{R}^3/\langle\sigma_a\rangle = \{p_3 \ge 0\}`.
    Adjoining the sphere's ideal :math:`p_1^2 + p_2^2 + p_3 = 1` and
    eliminating :math:`p_3 = 1 - p_1^2 - p_2^2` leaves the CLOSED UNIT
    DISK :math:`\{p_1^2 + p_2^2 \le 1\} = D^2`.

    ⚠ The dimension does NOT drop — :math:`\dim = 2 - 0 = 2`, because
    :math:`H` is finite.  That single fact is what makes this entry
    structurally unlike ``S^2/SO(2)`` (where :math:`2 - 1 = 1`) and is
    why it, and not the first entry, forced the chart-vs-section ruling:
    with no reduction, the chart buys nothing and the section is
    canonical.

    :math:`\det P = 4 x_a^2` vanishes exactly on the mirror's own
    fixed-point set, the great circle :math:`x_a = 0` — which in the
    realization's coordinates is the disk's BOUNDARY, and is a circle,
    not a finite point set.

    Full derivation with SymPy output, the Molien completeness check and
    the ``dim(m/m^2)`` minimality check:
    ``scratch/sigma_y_orbit_derivation.md`` (untracked; the load-bearing
    content is reproduced in ``docs/theory/foundations/manifolds.rst``).
    """
    import sympy as sp

    axis = group.mirror_axis
    if axis is None:
        # An admission contract, so a real raise: `-O` strips `assert`,
        # and this is the one precondition the derivation cannot check
        # from its own arguments.  Reachable only by registering this
        # builder under a non-Mirror key.
        raise ValueError(
            f"the reflection orbit-space derivation needs a mirror axis, "
            f"and {group.name!r} has none. Register "
            f"_sphere_mod_mirror only against Mirror entries."
        )
    coords = sp.symbols("x y z", real=True)
    x_a = coords[axis]
    kept = [c for i, c in enumerate(coords) if i != axis]
    u = sp.symbols("u0:3", real=True)

    invariants = (kept[0], kept[1], x_a**2)
    grad = [[sp.diff(q, v) for v in coords] for q in invariants]
    gram_xyz = sp.Matrix(
        [[sum(gi * gj for gi, gj in zip(a, b)) for b in grad] for a in grad]
    )
    # Re-express in the invariants: <grad p_3, grad p_3> = 4 x_a^2 = 4 p_3.
    gram = gram_xyz.subs({x_a**2: u[2]}).applyfunc(sp.simplify)

    return Quotient(
        base=base,
        by=group,
        # The CHART codomain, in invariant coordinates (p_1, p_2).
        realization=Ball(2),
        # The SECTION's image, in the base's ambient coordinates — what
        # DiscreteMeasure.quotient actually emits.  Closed, not strict:
        # [M] the cylindrical march seeds a level at x_a = 0 exactly.
        fundamental_domain=FundamentalDomain(
            base=base,
            normals=(tuple(1.0 if i == axis else 0.0 for i in range(3)),),
            label=f"{'xyz'[axis]}>=0",
        ),
        generators=(
            u[0] - invariants[0],
            u[1] - invariants[1],
            u[2] - invariants[2],
        ),
        syzygy=(),
        gram=gram,
        det_gram=sp.simplify(gram.det()),
        derived_by="hand",
        # det P = 4 p_3 = 4(1 - u0^2 - u1^2) after the sphere ideal: the
        # disk's boundary circle.
        singular_stratum=sp.simplify(1 - u[0] ** 2 - u[1] ** 2),
    )


def _mod_trivial(base: Manifold, group: "SubgroupOfO3") -> Quotient:
    r"""``M/{e} = M`` — the identity quotient, derived not tabulated.

    The invariant ring of the trivial group is all of
    :math:`\mathbb{R}[x_1..x_n]`, minimally generated by the coordinates,
    so :math:`P_{ij} = \langle \nabla x_i, \nabla x_j \rangle =
    \delta_{ij}` and :math:`\det P = 1`.  A determinant that vanishes
    nowhere means an empty singular stratum, i.e. a free action — which
    is right, vacuously, because the only group element is the identity.

    ⭐ This is the shipped ``AngularSymmetry.support``'s ``Trivial``
    answer, re-derived rather than tabulated: that property maps
    ``Trivial -> "S^2"`` in the string vocabulary and is the twin this
    type exists to absorb.
    """
    import sympy as sp

    n = _ambient(base)
    coords = sp.symbols(f"x0:{n}", real=True)
    return Quotient(
        base=base,
        by=group,
        realization=base,
        generators=tuple(coords),
        syzygy=(),
        gram=sp.eye(n),
        det_gram=sp.Integer(1),
        derived_by="hand",
        # M/{e} = M, so the identity map is a section and its image is all
        # of M — a fundamental domain cut by NO half-spaces.
        fundamental_domain=FundamentalDomain(base, (), "all"),
        # det P = 1 vanishes nowhere: no stratum, a free action (vacuously,
        # the only element being the identity).
        singular_stratum=None,
    )


#: ``(manifold type, group name) -> builder``.  Keyed on the PAIR because
#: a quotient is binary dispatch: it is a property of neither operand
#: alone.  About a dozen entries are expected; the engine that would
#: compute them instead of reading them is deferred, not refused.
_ORBIT_CATALOGUE: dict[tuple[type, str], Any] = {
    (Sphere, "SO2"): _sphere_mod_so2,
    # All three mirrors share ONE derivation — it reads the axis off the
    # group — so they are three keys, not three procedures.
    (Sphere, "sigma_x"): _sphere_mod_mirror,
    (Sphere, "sigma_y"): _sphere_mod_mirror,
    (Sphere, "sigma_z"): _sphere_mod_mirror,
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
        case RealSpace(d=d) | Ball(d=d):
            return d
        case FundamentalDomain(base=base):
            return _ambient(base)
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
