r"""Rigid motions of :math:`\mathbb{R}^d` — the Euclidean group :math:`E(d)`.

This module is the single place the project builds and composes geometric
transformations. Everything it exposes is ``numpy``-only: no transport
vocabulary, no measures, no neutrons. That is deliberate — a rotation is a
fact about space, not about what is being transported through it.

The algebra
-----------

The Euclidean group is a **semidirect** product,

.. math::

    E(d) \;=\; O(d) \ltimes \mathbb{R}^d ,

whose elements are pairs :math:`(Q, t)` acting on a point by
:math:`x \mapsto Qx + t`, with :math:`Q^\mathsf{T}Q = I` (hence
:math:`\det Q = \pm 1`). "Semidirect" rather than "direct" is the whole
content of the composition law: the left factor **acts on** the right one,

.. math::

    (Q_1, t_1)\circ(Q_2, t_2) = (Q_1 Q_2,\; Q_1 t_2 + t_1),
    \qquad
    (Q, t)^{-1} = (Q^\mathsf{T},\; -Q^\mathsf{T} t),

so translations do not simply add — the outer rotation carries the inner
translation with it. Getting this wrong is a whole bug class, which is why
composition is spelled once here (:meth:`RigidMotion.__matmul__`) and nowhere
else.

Why the affine part is not optional
-----------------------------------

A reflection or a rotation **in space** cannot be described without it. A
hyperplane is :math:`\{x : \hat n \cdot x = d\}`; only :math:`d = 0` passes
through the origin. A rotation needs a centre. Only the origin-fixing special
case is linear, and restricting to it is what made
``Mirror('x').is_invariant(mesh)`` read ``False`` for every production mesh —
the meshes start at :math:`x = 0` while their mirror plane sits at
:math:`x = a`, i.e. :math:`d \neq 0`, which was unexpressible.

The relation between the linear and affine cases is **one operation**:
conjugation by a translation (:meth:`RigidMotion.seated_at`),

.. math::

    Q \text{ seated at } c
    \;=\; T_c \circ Q \circ T_{-c}
    \;=\; \bigl(Q,\; (I - Q)c\bigr),

which reproduces every case — a reflection at signed offset :math:`d` is
:math:`(I - 2\hat n\hat n^\mathsf{T},\, 2d\,\hat n)`, a rotation about a centre
:math:`c` is :math:`(R,\,(I-R)c)`, an inversion through :math:`c` is
:math:`(-I,\,2c)`, and a translation is :math:`(I,\,v)`.

**Scope: rigid motions only** (:math:`Q^\mathsf{T}Q = I`), not general affine
maps — and this is a *separate* question from the one above, easy to run
together with it. The affine part is forced because you must be able to say
*where* a mirror is; excluding shear and scaling is a different claim, and it
does **not** follow from that proof. It rests on this:

**Orthogonality is not a convenience — it is what makes "symmetry" a
well-posed question.**

* A symmetry **preserves the structure**, and for a weighted point set that
  structure is distance and volume. A non-isometry carries the set to a
  *different* set, so :meth:`RigidMotion.permutes` and
  :meth:`RigidMotion.preserves` would not be answering a hard question — they
  would be answering a malformed one.
* A point group **is by definition** a subgroup of :math:`O(d)` seated at a
  point, not of :math:`GL(d)`.
* :func:`close_group` of a non-isometry is generically **infinite** — a scaling
  by 2 generates :math:`\mathbb{Z}` — so the finite-group machinery loses its
  meaning and :class:`NotAFinitePointGroupError` would fire on ordinary input.
* :meth:`RigidMotion.inverse` is :math:`Q^\mathsf{T}`, **exact**. A general
  linear part makes it a solve carrying conditioning, and every law that closes
  at round-off today would degrade to tolerance-bound.

Orthogonality is therefore the invariant enforced at construction, and a
non-isometry is not a value this type can hold. That the codebase happens to
have no consumer for shear is a footnote, never the argument.

**Where a genuinely non-rigid affine map would go — not here.** The real
candidate in this tree is the reference-cell :math:`\to` physical-cell map
behind the tensor-Legendre spatial basis, which carries an anisotropic scaling.
That is a **chart**: a different object with a different job, belonging beside
the curvilinear charts (:math:`(r,\theta) \leftrightarrow (x,y)`) in
:mod:`orpheus.geometry.coord`, which are **non-linear** and are likewise a
different concept rather than a degenerate case of this one.

Elements are parameterised by the complement of what they fix
-------------------------------------------------------------

============  ==========================  ==========================
              reflection                  rotation
============  ==========================  ==========================
fixed set     dimension :math:`d-1`       dimension :math:`d-2`
given by      its 1-D complement: a       its 2-D complement: a
              **normal**                  **rotation plane** + angle
:math:`\det`  :math:`-1`                  :math:`+1`
============  ==========================  ==========================

Both are therefore **dimension-generic**, and this module treats them so. In
:math:`\mathbb{R}^1` a reflection has normal :math:`[1]` and is
:math:`x \mapsto -x`; in :math:`\mathbb{R}^3` a rotation's plane has a normal
of its own, which is the familiar *axis* — a coincidence of :math:`d = 3`, not
a definition. Writing the rotation in its plane rather than about its axis is
what makes one formula serve :math:`d = 2` and :math:`d = 3` alike, and is why
:meth:`RigidMotion.rotation_about_axis` is a three-dimensional *convenience*
that delegates rather than a second implementation.

Points versus directions
------------------------

A point transforms affinely; a direction transforms by the linear part alone.
The two actions are therefore spelled separately (:meth:`RigidMotion.on_points`
and :meth:`RigidMotion.on_directions`) and there is no ``__call__`` that would
let a reader guess which is meant. The transport boundary law
:math:`\psi_{\text{in}}(x,\Omega) = \psi_{\text{out}}(g^{-1}x,\, Q_g^{-1}\Omega)`
reads directly in this vocabulary, and silently translating a direction — the
bug the split exists to prevent — cannot be written.
"""
from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # pragma: no cover - typing only
    from numpy.typing import ArrayLike

__all__ = [
    "NotAFinitePointGroupError",
    "Permutation",
    "RigidMotion",
    "close_group",
]


class NotAFinitePointGroupError(ValueError):
    r"""A generating set generates an INFINITE group.

    Named rather than bare, because the configuration that raises it is
    ordinary rather than exotic: **two parallel mirrors generate the infinite
    dihedral group**, and two parallel mirrors are the two reflective faces of
    a slab. Their product is the pure translation by twice their separation —
    an element of infinite order that a point group cannot contain.

    Every finite subgroup of :math:`E(d)` fixes a point (its own centroid), so
    this error is the type system saying "there is no common seat": either the
    generators are genuinely a space group rather than a point group, or one
    of them was built with an offset it should not carry.
    """


#: Tolerance on :math:`\|Q^\mathsf{T}Q - I\|_\infty` at construction.
#:
#: Directly-built elements (Householder, planar rotation, signed permutation)
#: are orthogonal to a few ULP; the slack is for accumulated products. Chosen
#: against measurement, not habit: closing :math:`I_h` (120 elements) reaches
#: `1.6e-15`, and a deliberately pathological 5000-deep chain of random
#: rotations reaches `2.8e-14` — so this leaves a factor of ~36 over a case no
#: real consumer produces, and ~640 over the deepest realistic one.
#:
#: It is a guard against a *genuine* non-isometry (a shear, a scaling), whose
#: departure is :math:`O(1)`. Note what that means for gates: an element with
#: a `1e-13` shear is a LEGAL value of this type, so a test asserting
#: orthogonality to `1e-14` on an arbitrary element would be asserting a
#: property the type does not promise. Gate the constructors' actual quality
#: separately from the type's invariant.
_ORTHOGONALITY_ATOL = 1e-12

#: Decimal places used to key group elements during closure.
#:
#: Distinct elements of the point groups this project closes are separated by
#: far more than this granularity — the closest pair in :math:`C_n` differs by
#: :math:`2\pi/n` — so rounding cannot merge two of them. A rounding boundary
#: could at worst *split* one element into two, which surfaces immediately as a
#: wrong group order rather than as a silent identification.
_CLOSURE_KEY_DECIMALS = 9

#: Refusal threshold for :func:`close_group`, so a generating set that does not
#: generate a finite group fails loudly instead of exhausting memory.
_MAX_GROUP_ORDER = 2000

#: Refusal threshold for :meth:`RigidMotion.element_order`.
_MAX_ELEMENT_ORDER = 1000

#: Singular-value threshold for the rank used by
#: :attr:`RigidMotion.fixed_subspace_dimension`.
_RANK_ATOL = 1e-10


@dataclass(frozen=True, eq=False)
class Permutation:
    r"""A bijection of :math:`\{0, \dots, n-1\}`, stored as the index array
    :math:`\pi` with :math:`i \mapsto \pi(i)`.

    Bijectivity is asserted at construction, so **returning one at all is the
    proof**. That is the point of the type: a nearest-neighbour match loop
    computes a *relation*, and a many-to-one relation satisfies every assertion
    a bare index array can carry — which is how ERR-073 shipped. Here the
    illegal state is unrepresentable rather than merely guarded against.

    Composition follows the maps it represents: :math:`(\pi \circ \rho)(i) =
    \pi(\rho(i))`, i.e. ``(a @ b).indices == a.indices[b.indices]``. Paired
    with :meth:`RigidMotion.permutes`' convention :math:`g(x_i) = x_{\pi(i)}`
    this makes :math:`\pi` a group **homomorphism**,
    :math:`\pi_{g \circ h} = \pi_g \circ \pi_h` — the one law that pins the
    composition order, the row-versus-column action, and the
    :math:`\pi`-versus-:math:`\pi^{-1}` choice all at once. None of the three
    is checkable alone.
    """

    indices: np.ndarray

    def __post_init__(self) -> None:
        pi = np.asarray(self.indices, dtype=np.int64)
        if pi.ndim != 1:
            raise ValueError(f"a permutation is a 1-D index array, got {pi.shape}")
        if not np.array_equal(np.sort(pi), np.arange(pi.size)):
            raise ValueError(
                f"indices are not a permutation of range({pi.size}): "
                f"{pi.size - np.unique(pi).size} value(s) repeat"
            )
        object.__setattr__(self, "indices", pi)

    @property
    def n(self) -> int:
        """The size of the set permuted."""
        return int(self.indices.size)

    @classmethod
    def identity(cls, n: int) -> "Permutation":
        return cls(np.arange(n, dtype=np.int64))

    def __matmul__(self, other: "Permutation") -> "Permutation":
        r""":math:`(\pi \circ \rho)(i) = \pi(\rho(i))` — ``other`` applied
        first."""
        if not isinstance(other, Permutation):
            return NotImplemented
        if other.n != self.n:
            raise ValueError(
                f"cannot compose permutations of {self.n} and {other.n} points"
            )
        return Permutation(self.indices[other.indices])

    def inverse(self) -> "Permutation":
        r""":math:`\pi^{-1}`, the unique bijection with
        :math:`\pi^{-1}\circ\pi = \mathrm{id}`."""
        out = np.empty_like(self.indices)
        out[self.indices] = np.arange(self.n, dtype=np.int64)
        return Permutation(out)

    @property
    def fixed_points(self) -> np.ndarray:
        r"""Indices with :math:`\pi(i) = i`.

        For a permutation induced by a group element this is exactly the set
        of points that element **stabilises**, so the union over the
        non-identity elements is the orbifold singular set :math:`\Sigma`.
        Membership is an INTEGER identity — no tolerance anywhere. The only
        place a tolerance ever enters is matching points while *building*
        :math:`\pi`, which is the one place the question is honestly numerical.
        """
        return np.flatnonzero(self.indices == np.arange(self.n))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Permutation):
            return NotImplemented
        return bool(np.array_equal(self.indices, other.indices))

    def __hash__(self) -> int:
        return hash(self.indices.tobytes())

    def __repr__(self) -> str:  # pragma: no cover - display only
        return f"Permutation({self.indices.tolist()})"


@dataclass(frozen=True, eq=False)
class RigidMotion:
    r"""An element :math:`(Q, t)` of :math:`E(d)`, acting as
    :math:`x \mapsto Qx + t`.

    Parameters
    ----------
    linear:
        The orthogonal part :math:`Q`, shape ``(d, d)``. Rejected at
        construction unless :math:`Q^\mathsf{T}Q = I` to
        :data:`_ORTHOGONALITY_ATOL` — the invariant lives here and only here,
        so every same-type-producing operation re-establishes it for free by
        routing back through construction.
    translation:
        The translation part :math:`t`, shape ``(d,)``. ``None`` means the
        origin-fixing (purely linear) element.

    Prefer the named constructors — :meth:`reflection`, :meth:`rotation`,
    :meth:`inversion`, :meth:`translation_by`, :meth:`signed_permutation` —
    which take the *geometric* data that determines the element rather than
    its matrix.
    """

    linear: np.ndarray
    translation: np.ndarray = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        Q = np.asarray(self.linear, dtype=float)
        if Q.ndim != 2 or Q.shape[0] != Q.shape[1]:
            raise ValueError(
                f"linear part must be a square matrix, got shape {Q.shape}"
            )
        d = Q.shape[0]
        departure = np.max(np.abs(Q.T @ Q - np.eye(d)))
        if departure > _ORTHOGONALITY_ATOL:
            raise ValueError(
                "linear part is not orthogonal: "
                f"max|QᵀQ − I| = {departure:.3e} > {_ORTHOGONALITY_ATOL:.1e}. "
                "RigidMotion models isometries only — a shear or a scaling is "
                "not a rigid motion and has no place in this type."
            )
        t = (
            np.zeros(d)
            if self.translation is None
            else np.asarray(self.translation, dtype=float)
        )
        if t.shape != (d,):
            raise ValueError(
                f"translation must have shape ({d},) to match the linear "
                f"part, got {t.shape}"
            )
        object.__setattr__(self, "linear", Q)
        object.__setattr__(self, "translation", t)

    # -- identity and structure ------------------------------------------

    @property
    def dimension(self) -> int:
        """The dimension :math:`d` of the space acted upon."""
        return int(self.linear.shape[0])

    @property
    def determinant(self) -> float:
        r""":math:`\det Q`, which is :math:`+1` (proper) or :math:`-1`
        (improper).

        This is *computed*, never declared. Whether an element is a rotation
        or a reflection is a property of the matrix, not a label attached to
        it — so there is no ``Reflection`` type here whose claim could drift
        from its content.
        """
        return float(np.linalg.det(self.linear))

    @property
    def is_proper(self) -> bool:
        r"""``True`` iff :math:`\det Q = +1` — orientation-preserving."""
        return self.determinant > 0.0

    @property
    def is_linear(self) -> bool:
        """``True`` iff the element fixes the origin (:math:`t = 0`)."""
        return bool(np.all(self.translation == 0.0))

    @property
    def is_translation(self) -> bool:
        r"""``True`` iff the linear part is the identity (:math:`Q = I`).

        The exact complement question to :attr:`is_linear` — that one asks
        whether the *translation* part is trivial, this one whether the
        *linear* part is. The identity motion answers ``True`` to both (the
        trivial translation), which is the group-theoretic reading: the
        translations are the kernel of :math:`(Q, t) \mapsto Q`, and the
        kernel contains the unit.

        Why this is the load-bearing predicate for the boundary tier: the
        action on **directions** is the linear part alone
        (:meth:`on_directions` drops :math:`t`), so "does this motion move
        angle at all?" IS this question. Both deck-element types derive
        their ``permutes_ordinates`` from it rather than declaring it —
        one source of truth for a fact the motion already knows.

        Exact compare, deliberately (like :attr:`is_linear`): the elements
        this discriminates are CONSTRUCTED (:meth:`identity`,
        :meth:`translation_by`, :meth:`reflection`, :meth:`rotation`), not
        measured, so a tolerance would only admit near-identities whose
        classification should be decided by the constructor that made them.
        """
        return bool(np.array_equal(self.linear, np.eye(self.dimension)))

    @property
    def fixed_subspace_dimension(self) -> int:
        r"""Dimension of :math:`\ker(Q - I)` — the linear part's fixed
        subspace.

        This is the structural claim that makes reflections and rotations
        dimension-generic: a reflection of :math:`\mathbb{R}^d` fixes a
        hyperplane (:math:`d-1`) and a simple rotation fixes a
        :math:`(d-2)`-dimensional subspace, for every :math:`d`.

        It reads the **linear** part, so it is unaffected by seating. A
        rotation by an angle small compared with :data:`_RANK_ATOL` is
        indistinguishable from the identity here, as it is numerically.
        """
        d = self.dimension
        rank = np.linalg.matrix_rank(self.linear - np.eye(d), tol=_RANK_ATOL)
        return d - int(rank)

    def element_order(self, *, atol: float = 1e-9) -> int | None:
        r"""The smallest :math:`k \ge 1` with :math:`g^k = e`, or ``None`` if
        the element has infinite order.

        Returns the order itself rather than answering "is the order
        :math:`n`?", because :math:`g^n = e` is satisfied by every element
        whose order *divides* :math:`n` — an involution passes a
        ":math:`g^4 = e`" check. Only the smallest such :math:`k` pins the
        element.

        ``None`` covers both irrational rotation angles and elements whose
        linear part has finite order but which carry a translation off the
        fixed subspace (a glide reflection, a screw): these generate infinite
        cyclic groups, and reporting the linear part's order for them would be
        a lie about the element.
        """
        identity = RigidMotion.identity(self.dimension)
        power = self
        for k in range(1, _MAX_ELEMENT_ORDER + 1):
            if power.approximately_equals(identity, atol=atol):
                return k
            power = self @ power
        return None

    # -- the group operations --------------------------------------------

    def __matmul__(self, other: "RigidMotion") -> "RigidMotion":
        r"""Composition :math:`(Q_1,t_1)\circ(Q_2,t_2)
        = (Q_1Q_2,\, Q_1 t_2 + t_1)`.

        ``a @ b`` is "apply ``b`` first, then ``a``", matching the reading of
        :math:`(a \circ b)(x) = a(b(x))` and of matrix products.
        """
        if not isinstance(other, RigidMotion):
            return NotImplemented
        if other.dimension != self.dimension:
            raise ValueError(
                "cannot compose rigid motions of different dimensions: "
                f"{self.dimension} and {other.dimension}"
            )
        return RigidMotion(
            self.linear @ other.linear,
            self.linear @ other.translation + self.translation,
        )

    def inverse(self) -> "RigidMotion":
        r""":math:`(Q, t)^{-1} = (Q^\mathsf{T},\, -Q^\mathsf{T}t)`.

        Uses :math:`Q^{-1} = Q^\mathsf{T}` directly. For an orthogonal matrix
        the transpose *is* the inverse exactly — no solve, no conditioning,
        and the result is bit-identical to the original's entries permuted.
        """
        Q_inv = self.linear.T
        return RigidMotion(Q_inv, -(Q_inv @ self.translation))

    def conjugated_by(self, other: "RigidMotion") -> "RigidMotion":
        r""":math:`h g h^{-1}` — the same transformation, seen in ``other``'s
        frame.

        Conjugation is how a transformation is *moved*: the conjugate of the
        reflection in the hyperplane :math:`H` by :math:`h` is the reflection
        in :math:`h(H)`. Both "re-seat this mirror somewhere else" and
        "re-orient it" are this one operation, which is why
        :meth:`seated_at` delegates here.
        """
        return other @ self @ other.inverse()

    def embedded_in(self, dimension: int) -> "RigidMotion":
        r"""The image of :math:`\iota : E(d) \hookrightarrow E(D)`,
        :math:`(Q, t) \mapsto (\mathrm{diag}(Q, I),\, (t, 0))`.

        The standard inclusion, and the operation that makes "dimension
        generic" a property **of the type** rather than a habit of its callers.
        Without it, relating a 1-D element to its 3-D counterpart is
        ``np.eye(3)`` padding at each call site — and then the 1-D/3-D split
        this module exists to delete has merely *moved*.

        It is a group homomorphism (:math:`\iota(g)\iota(h) = \iota(gh)`) that
        preserves the determinant and raises the fixed-subspace dimension by
        exactly :math:`D - d`. So a reflection stays a reflection and keeps
        fixing a hyperplane: in :math:`\mathbb{R}^1`, ``reflection(normal=[1])``
        fixes a 0-dimensional set (the origin); embedded in
        :math:`\mathbb{R}^3` it is :math:`\sigma_x`, fixing the plane
        :math:`x = 0` — the same element, read in a bigger room.

        The complementary axes are fixed pointwise, which is precisely the
        content of the :math:`(\mu, 0, 0)` embedding the angular layer already
        performs on its data.
        """
        d = self.dimension
        if dimension < d:
            raise ValueError(
                f"an embedding cannot lower dimension: {d} -> {dimension}. "
                "There is no canonical restriction; a projection would not be "
                "an isometry of the target."
            )
        Q = np.eye(dimension)
        Q[:d, :d] = self.linear
        t = np.zeros(dimension)
        t[:d] = self.translation
        return RigidMotion(Q, t)

    def seated_at(self, centre: "ArrayLike") -> "RigidMotion":
        r"""Conjugation by the translation to ``centre``:
        :math:`(Q,\, t + (I-Q)c)`.

        For a linear element this is exactly :math:`Q` *seated at* ``centre``
        — the same orthogonal transformation with ``centre`` in its fixed set
        instead of the origin. Applied to an element that already carries a
        translation it is the general conjugation, which is the honest
        generalisation and not idempotent.
        """
        return self.conjugated_by(RigidMotion.translation_by(centre))

    # -- the two actions --------------------------------------------------

    def on_points(self, points: "ArrayLike") -> np.ndarray:
        r"""The affine action :math:`x \mapsto Qx + t`.

        Accepts a single point of shape ``(d,)`` or a stack ``(..., d)`` and
        returns the same shape.
        """
        x = np.asarray(points, dtype=float)
        if x.shape[-1] != self.dimension:
            raise ValueError(
                f"points must have trailing dimension {self.dimension}, "
                f"got shape {x.shape}"
            )
        return x @ self.linear.T + self.translation

    def on_directions(self, directions: "ArrayLike") -> np.ndarray:
        r"""The linear action :math:`\Omega \mapsto Q\Omega`.

        A direction is not a point: it has no position, so the translation
        does not act on it. Applying the full affine map to a direction is a
        real bug class (it silently denormalises unit vectors); spelling the
        two actions separately makes it unwriteable.
        """
        v = np.asarray(directions, dtype=float)
        if v.shape[-1] != self.dimension:
            raise ValueError(
                f"directions must have trailing dimension {self.dimension}, "
                f"got shape {v.shape}"
            )
        return v @ self.linear.T

    def fixes(self, points: "ArrayLike", *, atol: float = 1e-12) -> np.ndarray:
        """Element-wise ``True`` where ``g(x) == x`` — the *affine* fixed set."""
        x = np.asarray(points, dtype=float)
        return np.all(np.abs(self.on_points(x) - x) <= atol, axis=-1)

    # -- action on a finite point set -------------------------------------

    def permutes(
        self, points: "ArrayLike", *, atol: float
    ) -> "Permutation | None":
        r"""The permutation :math:`\pi` with :math:`g(x_i) = x_{\pi(i)}`, or
        ``None`` if ``g`` does not map the set onto itself.

        Returns the permutation rather than a ``bool``, because a predicate
        that internally builds the permutation and throws it away *is* the
        missing primitive: every downstream reading — orbits, stabilisers, the
        singular set — is a reading of :math:`\pi`, and a returned permutation
        makes its own bijectivity assertable where a ``bool`` makes it
        unfalsifiable.

        Two guards, and the second is the one that is easy to omit:

        1. every image is within ``atol`` of some point of the set;
        2. the match is a **bijection**.

        Nearest-neighbour matching alone proves only that every image has
        *a* partner, which is strictly weaker than "the action permutes the
        set": two distinct sources may land on one target, leaving some point
        unmatched entirely. Appending a bit-identical duplicate of any point
        to an invariant set is the minimal witness — the duplicated position
        then carries twice the mass of its image, and every match map is
        non-injective, while an injectivity-free check certifies it.
        Since :math:`\pi` maps an :math:`n`-set to an :math:`n`-set,
        injectivity is equivalent to bijectivity, so counting distinct targets
        suffices. (ERR-073.)
        """
        x = np.asarray(points, dtype=float)
        if x.ndim != 2 or x.shape[1] != self.dimension:
            raise ValueError(
                f"points must have shape (n, {self.dimension}), "
                f"got {x.shape}"
            )
        n = x.shape[0]
        moved = self.on_points(x)
        # Nearest-neighbour match, all sources at once. Point counts here are
        # small (Lebedev-17 is 110), so the (n, n) distance matrix is cheaper
        # than a per-point Python loop — and one match rule cannot disagree
        # with itself the way a fast path plus a fallback can.
        dist = np.linalg.norm(moved[:, None, :] - x[None, :, :], axis=2)
        pi = np.argmin(dist, axis=1)
        if np.any(dist[np.arange(n), pi] > atol):
            return None  # some image is not a point of the set at all
        if np.unique(pi).size != n:
            return None  # the match is not a bijection (ERR-073)
        return Permutation(pi.astype(np.int64))

    def preserves(
        self,
        points: "ArrayLike",
        weights: "ArrayLike",
        *,
        atol: float,
        weight_atol: float,
    ) -> "Permutation | None":
        r"""The permutation :math:`\pi` with :math:`g(x_i) = x_{\pi(i)}`
        **and** :math:`w_i = w_{\pi(i)}`, or ``None``.

        :meth:`permutes` plus exactly one guard. The two are separate
        questions and so have separate names: a point set can be carried onto
        itself by ``g`` while the weights it carries are not — such a set is
        not invariant for any integration purpose, and asking one question
        with an optional ``weights=None`` would hide which was asked.

        The two tolerances are both required and both explicit because they
        are genuinely different windows. A coordinate is the accumulated
        result of a matrix product, so it carries more round-off than a
        weight, which is usually read straight from a table; a caller that
        wants them equal must say so.
        """
        pi = self.permutes(points, atol=atol)
        if pi is None:
            return None
        w = np.asarray(weights, dtype=float)
        if w.shape != (pi.n,):
            raise ValueError(
                f"weights must have shape ({pi.n},) to match the points, "
                f"got {w.shape}"
            )
        if np.any(np.abs(w[pi.indices] - w) > weight_atol):
            return None
        return pi

    # -- equality ---------------------------------------------------------

    def _exact_key(self) -> bytes:
        """Canonical bytes for exact equality and hashing.

        ``+ 0.0`` canonicalises ``-0.0`` to ``0.0`` so that the two compare
        AND hash alike; routing both dunders through this one function is
        what keeps that guarantee from drifting.
        """
        return (self.linear + 0.0).tobytes() + (self.translation + 0.0).tobytes()

    def __eq__(self, other: object) -> bool:
        """Exact (bit-level) equality. For the tolerant question, use
        :meth:`approximately_equals`."""
        if not isinstance(other, RigidMotion):
            return NotImplemented
        if other.dimension != self.dimension:
            return False
        return self._exact_key() == other._exact_key()

    def __hash__(self) -> int:
        return hash((self.dimension, self._exact_key()))

    def approximately_equals(
        self, other: "RigidMotion", *, atol: float = 1e-12
    ) -> bool:
        """``True`` iff both parts agree to ``atol``."""
        if other.dimension != self.dimension:
            return False
        return bool(
            np.all(np.abs(self.linear - other.linear) <= atol)
            and np.all(np.abs(self.translation - other.translation) <= atol)
        )

    def __repr__(self) -> str:  # pragma: no cover - display only
        kind = "proper" if self.is_proper else "improper"
        seated = "linear" if self.is_linear else f"t={self.translation}"
        return (
            f"RigidMotion(d={self.dimension}, {kind}, "
            f"fix={self.fixed_subspace_dimension}, {seated})"
        )

    # -- named constructors -----------------------------------------------

    @classmethod
    def identity(cls, dimension: int) -> "RigidMotion":
        """The identity of :math:`E(d)`."""
        return cls(np.eye(dimension))

    @classmethod
    def translation_by(cls, vector: "ArrayLike") -> "RigidMotion":
        r"""The pure translation :math:`(I, v)`.

        Named ``translation_by`` rather than ``translation`` because the
        latter is the field holding :math:`t`; a constructor and an accessor
        sharing a name is the kind of collision that survives review.
        """
        v = np.asarray(vector, dtype=float)
        if v.ndim != 1:
            raise ValueError(f"translation vector must be 1-D, got {v.shape}")
        return cls(np.eye(v.size), v)

    @classmethod
    def inversion(cls, dimension: int) -> "RigidMotion":
        r"""Inversion through the origin, :math:`x \mapsto -x`.

        Its determinant is :math:`(-1)^d`, so it is a reflection in
        :math:`\mathbb{R}^1`, a half-turn in :math:`\mathbb{R}^2`, and
        improper-but-not-a-reflection in :math:`\mathbb{R}^3` — one element,
        three familiar names, which is exactly why it is built from its
        definition rather than from any of them.
        """
        return cls(-np.eye(dimension))

    @classmethod
    def reflection(
        cls, *, normal: "ArrayLike", offset: float = 0.0
    ) -> "RigidMotion":
        r"""Reflection in the hyperplane
        :math:`\{x : \hat n \cdot x = \text{offset}\}`.

        .. math::

            x \;\mapsto\; x - 2(\hat n \cdot x - d)\,\hat n
            \;=\; (I - 2\hat n \hat n^\mathsf{T})\,x + 2d\,\hat n

        ``normal`` is keyword-only and is the **normal to the mirror**, never
        a direction lying in it. In :math:`\mathbb{R}^3` "reflection through
        the z axis" reads in English as :math:`\mathrm{diag}(-1,-1,1)`, which
        has :math:`\det = +1` and is a half-turn, not a reflection at all;
        forcing the caller to name the normal makes that reading unspellable
        rather than merely documented against.

        ``normal`` need not be a unit vector — it is normalised here — but a
        unit input is passed through exactly, so a coordinate normal yields a
        bit-exact signed-diagonal matrix.
        """
        n = np.asarray(normal, dtype=float)
        if n.ndim != 1:
            raise ValueError(f"normal must be a 1-D vector, got {n.shape}")
        length = float(np.linalg.norm(n))
        if length == 0.0:
            raise ValueError("normal must be non-zero — it names the mirror")
        n = n / length
        Q = np.eye(n.size) - 2.0 * np.outer(n, n)
        return cls(Q, 2.0 * float(offset) * n)

    @classmethod
    def rotation_from_circle_point(
        cls, *, plane: "ArrayLike", point: "ArrayLike"
    ) -> "RigidMotion":
        r"""Rotation in the 2-plane spanned by ``plane = (u, v)``, by the
        angle whose :math:`(\cos, \sin)` is ``point``.

        .. math::

            Q = I
              + \sin\theta\,(v u^\mathsf{T} - u v^\mathsf{T})
              + (\cos\theta - 1)(u u^\mathsf{T} + v v^\mathsf{T})

        acting as :math:`u \mapsto \cos\theta\,u + \sin\theta\,v`,
        :math:`v \mapsto -\sin\theta\,u + \cos\theta\,v`, and as the identity
        on the orthogonal complement — which is the :math:`(d-2)`-dimensional
        fixed set. One formula for every :math:`d \ge 2`.

        The primitive takes the **point on the unit circle**, not the angle,
        because the angle is a lossy parameterisation of it: the symmetric
        angles that matter here have exactly-representable cosines and sines
        (:math:`\cos(\pi/2) = 0` exactly) which ``np.cos`` of a rounded
        :math:`\pi/2` does not reproduce. A caller holding exact circle points
        — from a root-of-unity construction — keeps that exactness through
        this constructor. :meth:`rotation` is the convenience that computes
        the point from an angle.
        """
        P = np.asarray(plane, dtype=float)
        if P.shape[0] != 2 or P.ndim != 2:
            raise ValueError(
                f"plane must be two vectors, shape (2, d), got {P.shape}"
            )
        if P.shape[1] < 2:
            raise ValueError(
                "there is no rotation of R^1: SO(1) = {e}, so the only "
                "1-dimensional rotation is the identity — use "
                "RigidMotion.identity(1). A rotation fixes a (d-2)-dimensional "
                "subspace, and d-2 = -1 here. R^1's non-trivial content lives "
                "entirely in reflection()."
            )
        u, v = P[0], P[1]
        gram = P @ P.T
        if not np.allclose(gram, np.eye(2), atol=_ORTHOGONALITY_ATOL):
            raise ValueError(
                "the rotation plane must be given by an ORTHONORMAL pair; "
                f"max|GG^T − I| = {np.max(np.abs(gram - np.eye(2))):.3e}"
            )
        c, s = np.asarray(point, dtype=float)
        radius = float(np.hypot(c, s))
        if abs(radius - 1.0) > _ORTHOGONALITY_ATOL:
            raise ValueError(
                f"point must lie on the unit circle, got |({c}, {s})| = "
                f"{radius}"
            )
        Q = (
            np.eye(u.size)
            + s * (np.outer(v, u) - np.outer(u, v))
            + (c - 1.0) * (np.outer(u, u) + np.outer(v, v))
        )
        return cls(Q)

    @classmethod
    def rotation(cls, *, plane: "ArrayLike", angle: float) -> "RigidMotion":
        r"""Rotation in the 2-plane ``plane = (u, v)`` by ``angle``, taking
        :math:`u` toward :math:`v`.

        Delegates to :meth:`rotation_from_circle_point`; see there for the
        formula and for why the circle point is the primitive.
        """
        return cls.rotation_from_circle_point(
            plane=plane, point=(np.cos(angle), np.sin(angle))
        )

    @classmethod
    def rotation_about_axis(
        cls, *, axis: "ArrayLike", angle: float
    ) -> "RigidMotion":
        r"""Rotation of :math:`\mathbb{R}^3` about ``axis`` by ``angle``,
        right-handed.

        A three-dimensional convenience, not a second construction: it builds
        the rotation plane orthogonal to ``axis`` and delegates to
        :meth:`rotation`. "About an axis" is available only because in
        :math:`d = 3` the :math:`(d-2)`-dimensional fixed set is a line and a
        2-plane is named by its normal; in every other dimension the plane is
        the only spelling.

        The in-plane basis is chosen deterministically (seeded from the
        coordinate axis least aligned with ``axis``), so the returned matrix
        is reproducible. Any other admissible basis gives the same rotation
        to round-off, since a rotation about an axis by an angle is unique.
        """
        a = np.asarray(axis, dtype=float)
        if a.shape != (3,):
            raise ValueError(
                f"rotation_about_axis is three-dimensional; got axis shape "
                f"{a.shape}. In other dimensions name the rotation PLANE."
            )
        length = float(np.linalg.norm(a))
        if length == 0.0:
            raise ValueError("axis must be non-zero — it names the rotation")
        a = a / length
        seed = np.eye(3)[int(np.argmin(np.abs(a)))]
        u = seed - (seed @ a) * a
        u = u / np.linalg.norm(u)
        v = np.cross(a, u)
        return cls.rotation(plane=(u, v), angle=angle)

    @classmethod
    def signed_permutation(
        cls, *, permutation: Sequence[int], signs: Sequence[float]
    ) -> "RigidMotion":
        r"""The element sending :math:`e_{\text{permutation}[i]}` to
        :math:`\text{signs}[i]\,e_i`.

        The hyperoctahedral group :math:`(\mathbb{Z}_2)^d \rtimes S_d` — the
        symmetry group of the :math:`d`-cube, which is :math:`O_h` at
        :math:`d = 3`. Its elements have exactly one non-zero entry per row
        and column, both :math:`\pm 1`, so they are exactly orthogonal with no
        round-off at all, and they are the natural spelling for a cubic point
        group: a rotation-and-reflection decomposition of the same 48 matrices
        would introduce transcendental entries for elements whose true entries
        are integers.
        """
        perm = tuple(int(p) for p in permutation)
        d = len(perm)
        if sorted(perm) != list(range(d)):
            raise ValueError(
                f"permutation must be a permutation of range({d}), got {perm}"
            )
        sgn = np.asarray(signs, dtype=float)
        if sgn.shape != (d,) or not np.all(np.abs(sgn) == 1.0):
            raise ValueError(
                f"signs must be {d} values, each exactly ±1, got {signs!r}"
            )
        Q = np.zeros((d, d))
        for i, p in enumerate(perm):
            Q[i, p] = sgn[i]
        return cls(Q)


def close_group(
    generators: Iterable[RigidMotion],
    *,
    dimension: int | None = None,
) -> tuple[RigidMotion, ...]:
    r"""The group generated by ``generators`` — closure under composition.

    A realization is normally a *generating* set, not the group: the standard
    :math:`D_{nh}` generators number :math:`2n+1` for a group of order
    :math:`4n`. Checking that a set is carried onto itself needs only
    generators, but asking whether one group **contains** another needs the
    elements, so the closure is taken here rather than at each such question.

    Raises :class:`NotAFinitePointGroupError` if the generated group exceeds
    :data:`_MAX_GROUP_ORDER`, which is what an infinite group looks like from
    inside a breadth-first closure. That is not an exotic input: two PARALLEL
    mirrors compose to the pure translation by twice their separation, and two
    parallel mirrors are the two reflective faces of a slab.

    ``dimension`` is required only when ``generators`` is empty, where it
    fixes which trivial group is meant.
    """
    gens = list(generators)
    if not gens:
        if dimension is None:
            raise ValueError(
                "close_group of an empty generating set needs an explicit "
                "dimension — it cannot say which trivial group is meant"
            )
        return (RigidMotion.identity(dimension),)

    d = gens[0].dimension
    if dimension is not None and dimension != d:
        raise ValueError(
            f"dimension={dimension} contradicts the generators' dimension {d}"
        )
    if any(g.dimension != d for g in gens):
        raise ValueError("all generators must act on the same dimension")

    def key(g: RigidMotion) -> bytes:
        # Membership is a dict lookup on rounded entries, not a scan over the
        # accumulated elements: the closure costs O(|G|²) products, and an
        # O(|G|) scan inside each would make it cubic.
        return (
            (np.round(g.linear, _CLOSURE_KEY_DECIMALS) + 0.0).tobytes()
            + (np.round(g.translation, _CLOSURE_KEY_DECIMALS) + 0.0).tobytes()
        )

    identity = RigidMotion.identity(d)
    elements: list[RigidMotion] = [identity]
    seen: set[bytes] = {key(identity)}

    for g in gens:
        if key(g) not in seen:
            elements.append(g)
            seen.add(key(g))

    frontier = list(elements)
    while frontier:
        fresh: list[RigidMotion] = []
        for a in frontier:
            for b in list(elements):
                c = a @ b
                k = key(c)
                if k not in seen:
                    elements.append(c)
                    seen.add(k)
                    fresh.append(c)
                    if len(elements) > _MAX_GROUP_ORDER:
                        raise NotAFinitePointGroupError(
                            "generating set does not generate a finite group "
                            f"(exceeded {_MAX_GROUP_ORDER} elements). Every "
                            "finite subgroup of E(d) fixes a point, so this "
                            "generating set has no common seat — check for an "
                            "element carrying a translation off its own fixed "
                            "subspace (two parallel mirrors are the common "
                            "case: their product is a pure translation)."
                        )
        frontier = fresh
    return tuple(elements)
