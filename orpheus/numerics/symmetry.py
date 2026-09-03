r"""Subgroups of :math:`O(3)` — named, realized, and ordered by containment.

Quadrature *selection* is fundamentally a symmetry-matching problem.
A geometry has a natural symmetry group :math:`G_{\text{geom}}`
(slab :math:`\to` :math:`O(2)` about the slab normal, sphere :math:`\to`
:math:`O(3)`, hexagonal lattice :math:`\to` :math:`D_{6h}`, …) and an
angular quadrature carries an *invariance group*
:math:`G_{\text{quad}}` — the maximal subgroup of :math:`O(3)` that
permutes the quadrature nodes among themselves while preserving
weights. The selection rule is

.. math::
   :label: subgroup-of-o3-containment

   G_{\text{geom}} \subseteq G_{\text{quad}}
   \;\Longleftrightarrow\;
   \forall\, g \in G_{\text{geom}},\; g \in G_{\text{quad}}.

The "if and only if" is a tautology of subgroup containment, but its
*application* is non-trivial: given a geometry-quadrature pair, decide
whether the quadrature respects every symmetry the geometry exhibits.
A failure (e.g. integrating an :math:`O_h`-symmetric flux against an
:math:`SO(2)`-only product quadrature) introduces spurious
azimuthal-symmetry breaking that pollutes higher-order moments — the
classical observation behind Lebedev grids in computational chemistry
(Lebedev 1976, Stiefel & Fässler 1979).

ORPHEUS's relevant subgroup lattice is **finite and small** — six named
groups plus the parameterised families :math:`C_n`, :math:`D_{nh}`,
:math:`\sigma_a`, :math:`SO(2)_a` and :math:`O(2)_a` — and every relation
in it is COMPUTED, never tabulated. :math:`\mathfrak{so}(3)` is simple and
three-dimensional, so its Lie subalgebras are :math:`\{0\}`, one line
:math:`\mathbb{R}\,[\hat a]_\times` per axis :math:`\hat a`, and
:math:`\mathfrak{so}(3)` itself — there is no two-dimensional one. A closed
subgroup of :math:`O(3)` is therefore exactly the pair (identity component,
one representative per connected component), :class:`Realization`, and
containment, the normaliser, the identity component and the question
"does the connected part fix these nodes" are each ONE computation on
that pair (:class:`IdentityComponent`). The tag — :class:`_NamedSubgroup`,
:class:`Cn`, :class:`Dnh`, :class:`Mirror`, :class:`SO2`, :class:`O2` — is
the group's IDENTITY and name, and :func:`_realize` is the one place it is
read to build the realization. A hand-written relation TABLE and a
hand-written dispatch ARM lived here until 2026-09-02, and each shipped a
false edge before it went (``_NAMED_LATTICE``'s ``D_nh ⊆ O2``, retired
2026-08-02; ``_axial_contains``'s ``SO2('x') ⊉ C_1``, retired 2026-09-02):
the argument for computing rather than hand-maintaining, made twice, in
two different shapes.

Whether a MEASURE is invariant under a group is the measure's question, and
it lives with the measure: :meth:`DiscreteMeasure.is_invariant_under
<orpheus.numerics.measure.DiscreteMeasure.is_invariant_under>` and its siblings
delegate to :mod:`orpheus.numerics.invariance` (since R2 of #434, 2026-09-03 —
until then the kernel, the closure and the lattice walk lived here, this
module imported :mod:`~orpheus.numerics.measure` at module scope, and the
manifold module could reference a group only under ``TYPE_CHECKING``).  This
module imports :mod:`orpheus.geometry.transformation` and, inside
``numerics``, only :mod:`~orpheus.numerics.roots_of_unity` — itself a
numpy-only leaf (`[M]` 2026-09-03: its imports are ``numpy`` and
``numpy.typing``), and it must STAY one, because the day it imports
``measure`` or ``manifold`` this module re-acquires the cycle R2 broke: a
group is a group.

References
----------

* Hamermesh, M. (1962). *Group Theory and Its Application to Physical
  Problems*. Addison-Wesley. §2.5 (finite point groups), §9.4
  (octahedral and icosahedral groups).
* Stiefel, E. and Fässler, A. (1979). *Group Theoretical Methods and
  Their Applications*. Birkhäuser. Chapters 4-5 (crystallographic
  point groups, invariant theory of finite groups).
* Lebedev, V.I. (1976). "Quadratures on a sphere." *USSR Computational
  Mathematics and Mathematical Physics* **16**(2), 10-24. The
  octahedral-invariant construction the :math:`O_h` invariance check
  (:mod:`~orpheus.numerics.invariance`) validates.
* Carlson, B.G. and Lathrop, K.D. (1968). "Transport theory: the
  method of discrete ordinates." In *Computing Methods in Reactor
  Physics*, Greenspan, Kelber, Okrent, eds., Gordon & Breach.
  Level-symmetric :math:`S_N` quadratures, also :math:`O_h`-invariant
  by construction.

See Also
--------

:ref:`discrete-measures` (theory page) — section "Symmetry groups for
quadrature invariance" — narrates how invariance ties into the
quadrature-selection logic, which ships:
:func:`~orpheus.numerics.quadrature.registry.select_quadrature` consults
:meth:`SubgroupOfO3.contains` at stage 0 and the measure's own nodes at
stage 1.  (This paragraph said "Issue 5 … will install" until 2026-09-03.)
"""

from __future__ import annotations

import functools
from dataclasses import dataclass
from enum import Enum
from typing import ClassVar

import numpy as np

from orpheus.geometry.transformation import (
    RigidMotion,
    close_group as _close_rigid_motions,
)

from .roots_of_unity import roots_of_unity

#: The ambient coordinate an axis label names.  :math:`\mathbb{R}^3`'s axes
#: are the convention every realization below is built in (principal axis
#: along z, vertex on x), and the manifold module reads it — the direction
#: the runtime import graph runs since R2 of #434 (``manifold -> symmetry``;
#: until then it ran the other way and the table lived in ``manifold``, and
#: `[M]` re-importing it from there under the new direction closed a 2-cycle
#: that killed 6 of 9 entry points).  The mirror family names its plane by
#: the NORMAL, so ``AXIS_INDEX["x"]`` selects the normal of :math:`\sigma_x`
#: as well as the axis of :math:`SO(2)_x` and :math:`O(2)_x`.
AXIS_LETTER: tuple[str, str, str] = ("x", "y", "z")
#: The inverse, ``AXIS_INDEX[AXIS_LETTER[i]] == i`` — one spelling of the
#: convention, DERIVED from the letters so the two cannot disagree.
AXIS_INDEX: dict[str, int] = {a: i for i, a in enumerate(AXIS_LETTER)}

#: The absolute tolerance under which two realized operators, or a vector and
#: its image, are the same ELEMENT — the one band for every element-level
#: comparison a :class:`Realization` makes (membership, containment, the
#: normaliser).  Absolute, never relative: `[M]` ``np.allclose``'s default
#: ``rtol=1e-5`` on a unit component turned a written ``atol=1e-12`` into an
#: effective 1e-5 (elegance review, 2026-09-02).  The realized worst near-fix
#: residual over every shipped finite tag is 6.7e-16 and the realized
#: elements' orthogonality defect 1.3e-15 (qa census, 2026-09-02), so the band
#: is slack by six orders either way; what matters is that it is ONE band.
#: (Named ``_MEMBERSHIP_ATOL`` until 2026-09-02, when the same name in
#: ``manifold.py`` — a POINT on a manifold, 1e-12 — made a grep answer twice.)
_ELEMENT_ATOL: float = 1e-9

#: :math:`I_3`, minted once — every zero-dimensional membership test compares to it.
_IDENTITY_3: np.ndarray = np.eye(3)


# ---------------------------------------------------------------------------
# Group enumeration
# ---------------------------------------------------------------------------


class _NamedSubgroup(Enum):
    """Named, parameter-free subgroups of :math:`O(3)`.

    Parameterised families (:class:`Cn`, :class:`Dnh`, :class:`Mirror`,
    :class:`SO2`, :class:`O2`) are separate types so their parameter — an
    order, or an axis — is type-safe without dragging parameters into the
    enum machinery. Two entries have LEFT this enum for that reason, each
    after shipping a wrong answer: ``Z2`` became ``Mirror(axis)`` on
    2026-08-02 and ``SO2`` became ``SO2(axis)`` on 2026-09-01 (tracker
    2.4 of the angular-spaces campaign, #429).
    """

    Trivial = "Trivial"        # {e}
    Dinfh = "Dinfh"            # D_∞h — the FULL cylindrical group (see below)
    OctahedralOh = "Oh"        # full octahedral group, 48 elements
    IcosahedralIh = "Ih"       # full icosahedral group, 120 elements
    SO3 = "SO3"                # proper rotations of the sphere
    O3 = "O3"                  # orthogonal group of R^3 (top of lattice)

    # ``Dinfh`` was named ``O2`` until 2026-08-02, and the name was
    # wrong in a load-bearing way. The entry is REALIZED as axial
    # rotations + σ_h, which is C_∞h; true O(2) embedded in 3-D is
    # C_∞v (rotations + VERTICAL mirrors) — and since 2026-09-02 THAT
    # group ships as the axis-parameterised :class:`O2` (#432): the
    # stabiliser of the axis vector, which is what the old name should
    # have meant. Neither D_∞h nor O(2)_a contains D_{nh}, because
    # D_{nh} carries C₂ axes lying IN the plane, and those flip the
    # axis — so the lattice's asserted ``D_nh ⊆ O2`` was false under
    # both readings, and a committed test pinned it.
    #
    # D_∞h — C_∞ rotations, the in-plane C₂ axes, σ_h, the vertical
    # mirrors and the improper products — DOES contain every D_{nh},
    # so the relation the test asserted is correct once the name is.
    # (It is the continuous symmetry an axisymmetric cylinder HAS; the
    # geometry table does not spend it — `[M]` 2026-09-02 the cylinder
    # row's ``continuous_isotropy`` is ``Trivial`` and nothing outside
    # this module and its tests consumes ``Dinfh``.) In the new
    # vocabulary it is the stabiliser of the z-axis as a LINE,
    # D_∞h = O(2)_z × {e, σ_h}, and ``Dinfh.contains(O2("z"))`` is the
    # edge that says so.


@dataclass(frozen=True)
class Cn:
    """Cyclic group :math:`C_n` of :math:`n`-fold proper rotations
    about a fixed axis.

    Parameters
    ----------
    n : int
        Order of the cyclic group; must be :math:`\\ge 1`. ``Cn(1)``
        is the trivial group; ``Cn(2)`` is the rotation-only sibling
        of :class:`Mirror` — ``Cn(2)`` is specifically the
        :math:`180^{\\circ}`-ROTATION group (:math:`\\det = +1`), where a
        mirror is a reflection (:math:`\\det = -1`). Same order, different
        group, and the two spellings exist to keep them apart.
    """

    n: int

    def __post_init__(self) -> None:  # pragma: no cover - guard
        if self.n < 1:
            raise ValueError(f"Cn requires n >= 1, got n={self.n}")


@dataclass(frozen=True)
class Dnh:
    """Dihedral group :math:`D_{nh}` of order :math:`4n`.

    :math:`D_{nh}` is the full symmetry group of an :math:`n`-prism:
    :math:`n` proper rotations about the principal axis, :math:`n`
    rotations through perpendicular :math:`C_2` axes, plus a horizontal
    mirror plane and the improper rotations that combine these.

    Parameters
    ----------
    n : int
        Order of the principal axis; must be :math:`\\ge 1`. Common
        cases for ORPHEUS: ``Dnh(2)`` for rectangular Cartesian 2-D
        geometries, ``Dnh(4)`` for square lattices, ``Dnh(6)`` for
        hexagonal lattices (forward-looking — not yet consumed).
    """

    n: int

    def __post_init__(self) -> None:  # pragma: no cover - guard
        if self.n < 1:
            raise ValueError(f"Dnh requires n >= 1, got n={self.n}")


@dataclass(frozen=True)
class Mirror:
    r"""The order-2 reflection group :math:`\{e, \sigma_a\}`.

    **Parameterised by the plane, because the plane is not a
    convention.** This family replaced a parameter-free ``Z2`` entry on
    2026-08-02, and the entry was wrong in a way no test could see.

    ``Z2`` was realized as :math:`\sigma_z` while its docstring said
    *"pick* :math:`z \to -z` *as a canonical representative; any single
    reflection works, the choice is convention"*. `[M]` That is false on
    a rule this tree ships: ``product(4, 3)`` is closed under
    :math:`\sigma_z` and **not** under :math:`\sigma_x`, and ``Z2``
    answered ``True``. The claim survived because every fixture that
    reached it was symmetric in all three planes, where the answer is
    ``True`` either way — a degenerate fixture masking an overload, one
    more instance of the pattern this campaign keeps finding.

    Worse, the tag meant two DIFFERENT things by node shape: the 3-D arm
    tested :math:`\sigma_z`, the 1-D arm tested :math:`x \to -x` on the
    single coordinate (plane-free). The tree's canonical embedding of a
    polar marginal is :math:`(\mu, 0, 0)` (see
    :meth:`~orpheus.numerics.quadrature.directional.Quadrature.axis_cosines`),
    so the slab mirror :math:`\mu \to -\mu` **is** :math:`\sigma_x` — a
    different matrix from the one the 3-D arm used. `[M]` On an
    ASYMMETRIC :math:`\mu`-set the two arms disagree: the 3-D arm
    certifies ``True`` for a set that violates :math:`\mu \to -\mu`,
    which is a false certification in the dangerous direction (the
    ERR-072 family).

    Naming the plane is what makes the two arms say the same thing. On a
    1-D measure embedded as :math:`(\mu, 0, 0)`, :math:`\sigma_y` and
    :math:`\sigma_z` fix every node **pointwise** — so they hold
    trivially, and :math:`\sigma_x` is the one real test. That is now
    derived from the embedding rather than asserted by a branch.

    Parameters
    ----------
    axis : str
        The mirror's NORMAL, one of ``"x"`` / ``"y"`` / ``"z"``; the
        plane is the one perpendicular to it through the origin. So
        ``Mirror("z")`` is :math:`\sigma_h` for an axial geometry and
        ``Mirror("x")`` is the slab's polar mirror.

    See Also
    --------
    Cn
        The PROPER order-2 sibling is ``Cn(2)`` — a rotation, not a
        reflection. Keeping the two spellings distinct is the point:
        ``Mirror`` has :math:`\det = -1` and therefore lies in no
        proper-rotation group.
    """

    axis: str

    def __post_init__(self) -> None:  # pragma: no cover - guard
        if self.axis not in ("x", "y", "z"):
            raise ValueError(
                f"Mirror requires axis in x/y/z, got {self.axis!r}. The "
                f"axis names the mirror's NORMAL; there is no unnamed "
                f"reflection, which is exactly what the retired Z2 tag "
                f"pretended there was."
            )


@dataclass(frozen=True)
class SO2:
    r"""The continuous group :math:`C_\infty` of proper rotations about a
    named coordinate axis.

    **Parameterised by the axis, because the axis is not a convention** —
    the ruling :class:`Mirror` records for the plane, reached the same
    way, one month later. Until 2026-09-01 this was the parameter-free
    named entry ``SO2``, realized about :math:`z`: its exactness test
    asked whether every node had :math:`\rho = \sqrt{x^2 + y^2} = 0`.
    The slab's polar marginal embeds along :math:`x`, so `[M]`
    ``SO2.is_invariant`` read ``False`` on the one shipped rule whose
    orbit space IS :math:`S^2/SO(2)` — the axis-convention collision
    the angular-spaces plan (#429) recorded as its Part IV obstacle 1.

    The tree genuinely carries two poles, so no fiat can settle it. The
    slab and sphere geometries, and the real spherical-harmonic basis
    (``cos θ = μ_x`` in ``_evaluate_real_sh``), are about :math:`x`;
    the product rules, :math:`C_n`, :math:`D_{nh}` and
    :math:`D_{\infty h}` are about :math:`z`. One Gauss–Legendre rule on
    :math:`\mu` serves BOTH — as the slab's rule on :math:`S^2/SO(2)_x`
    and as the polar factor of a product on :math:`S^2/SO(2)_z` — so the
    group a marginal was quotiented by cannot be spelled without its axis,
    and two orbit spaces that differ only in the axis must not compare
    equal.

    Parameters
    ----------
    axis : str
        The rotation axis, one of ``"x"`` / ``"y"`` / ``"z"``.
        ``SO2("x")`` is the group the slab and sphere geometries spend
        (:data:`~orpheus.numerics.quadrature.registry.GEOMETRY_ANGULAR_SYMMETRY`).
        ``SO2("z")`` is the one the finite families sit inside —
        :math:`C_n \subset SO(2)_z \subset D_{\infty h}` — because those
        are realized about :math:`z` in the standard setting; a
        :math:`C_n` lies in no other axis's rotation group.

    See Also
    --------
    Mirror
        The reflection family, parameterised for the same reason.
    Cn
        The finite cyclic subgroups, about :math:`z` only. They are not
        axis-parameterised: no consumer has yet needed a :math:`C_n`
        about another axis, and the project unifies after the second
        instance, not before it.
    """

    axis: str

    def __post_init__(self) -> None:  # pragma: no cover - guard
        if self.axis not in ("x", "y", "z"):
            raise ValueError(
                f"SO2 requires axis in x/y/z, got {self.axis!r}. The axis "
                f"names the rotation axis; there is no unnamed axial "
                f"rotation group, which is exactly what the retired bare "
                f"SO2 entry pretended there was."
            )


@dataclass(frozen=True)
class O2:
    r"""The continuous group :math:`O(2)_a = C_{\infty v}` — the stabiliser
    of a named coordinate axis in :math:`O(3)`.

    Every rotation about the axis :math:`a` **and** every reflection in a
    plane CONTAINING it:

    .. math:: O(2)_a \;=\; \{\, g \in O(3) : g\,\hat e_a = \hat e_a \,\},

    the subgroup fixing the axis vector pointwise. Its identity component
    is :class:`SO2` (the proper half); the other component is the coset
    :math:`SO(2)_a\,\sigma_v` of the vertical mirrors. Embedded in three
    dimensions it is the point group :math:`C_{\infty v}` — NOT
    :math:`C_{\infty h}` (rotations plus the mirror PERPENDICULAR to the
    axis, which flips :math:`\hat e_a`), and not :math:`D_{\infty h}`,
    which is the stabiliser of the axis as a LINE,
    :math:`\{g : g\hat e_a = \pm\hat e_a\} = O(2)_a \times \{e, \sigma_h\}`.

    **Why it exists (#432, 2026-09-02): an orbit space is named by the
    LARGEST group with its orbits.** The :math:`SO(2)_a`-orbits on
    :math:`S^2` are the circles of constant :math:`\mu = \Omega \cdot
    \hat e_a`, and the group carrying every one of those circles onto
    itself is exactly the stabiliser of :math:`\hat e_a` — a vertical
    mirror reflects each circle into itself. So the invariant rings
    coincide,

    .. math:: \mathbb{R}[x]^{SO(2)_a} \;=\; \mathbb{R}[x]^{O(2)_a}
              \;=\; \mathbb{R}[\,x_a,\; x_b^2 + x_c^2\,],

    the orbit-space derivation is one procedure with one output, and the
    catalogue entry :math:`S^2/O(2)_a` records THIS group
    (``manifold._sphere_mod_o2``; asking for the quotient by
    :math:`SO(2)_a` is refused with the theorem, at the catalogue door and
    at :class:`~orpheus.numerics.manifold.Quotient`'s construction, both
    reading :attr:`SubgroupOfO3.orbit_stabiliser`). A basis on that
    entry — the Legendre polynomials :math:`P_\ell(\mu)` — therefore
    declares the FULL group its functions have, and the frame's lattice
    gate admits it on a :math:`\sigma_b`-folded rule (:math:`b \ne a`)
    because :math:`\sigma_b \in O(2)_a`. Until this member existed
    :attr:`~orpheus.numerics.basis.base.Basis.invariance_group` could only
    be DERIVED as the strict lower bound :math:`SO(2)_a`, and that
    mathematically admissible pairing was over-refused.

    ⚠ **The name has history.** ``O2`` was the name of the
    :math:`D_{\infty h}` entry until 2026-08-02 and was retired because
    that entry realizes :math:`C_{\infty h}` (see the note on
    :class:`_NamedSubgroup`). This class is the *other* group — the one the
    old name should have meant — and it is axis-parameterised for the
    reason :class:`SO2` is: the tree carries two poles.

    **Exact invariance.** A finite point set is :math:`O(2)_a`-closed iff
    every node lies ON the axis — the same criterion as :math:`SO(2)_a`,
    since :math:`SO(2)_a \subseteq O(2)_a` forces axis support and a point
    on the axis is fixed by every vertical mirror — decided, never sampled
    (ERR-072).

    Parameters
    ----------
    axis : str
        The fixed axis, one of ``"x"`` / ``"y"`` / ``"z"``. ``O2("x")`` is
        the group the slab and sphere geometries spend
        (:data:`~orpheus.numerics.quadrature.registry.GEOMETRY_ANGULAR_SYMMETRY`):
        a slab is invariant under every rotation about its normal AND
        under :math:`y \to -y`, which is why its flux depends on
        :math:`\Omega` only through :math:`\mu`.

    See Also
    --------
    SO2
        The proper half — the identity component.
    """

    axis: str

    def __post_init__(self) -> None:  # pragma: no cover - guard
        if self.axis not in ("x", "y", "z"):
            raise ValueError(
                f"O2 requires axis in x/y/z, got {self.axis!r}. The axis "
                f"names the axis the group fixes; there is no unnamed "
                f"axial group, for the reason there is no unnamed axial "
                f"rotation group."
            )


# ---------------------------------------------------------------------------
# Public façade
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class SubgroupOfO3:
    """Subgroup-of-:math:`O(3)` lattice with containment + invariance checks.

    Class attributes (named entries) are pre-instantiated singletons:

    .. code-block:: python

        SubgroupOfO3.Trivial
        SubgroupOfO3.Dinfh
        SubgroupOfO3.OctahedralOh
        SubgroupOfO3.IcosahedralIh
        SubgroupOfO3.SO3
        SubgroupOfO3.O3

    Parameterised families are constructors:

    .. code-block:: python

        SubgroupOfO3.Cn(6)        # cyclic order 6
        SubgroupOfO3.Dnh(6)       # dihedral D_6h (hex lattice)
        SubgroupOfO3.Mirror('z')  # {e, sigma_z} — the plane is named
        SubgroupOfO3.SO2('x')     # rotations about x — the axis is named
        SubgroupOfO3.O2('x')      # rotations about x AND the mirrors through it

    Containment via :meth:`contains` implements
    :eq:`subgroup-of-o3-containment`, COMPUTED from the two groups'
    :attr:`realization`.  Whether a measure is invariant under this group is
    the MEASURE's question — :meth:`DiscreteMeasure.is_invariant_under
    <orpheus.numerics.measure.DiscreteMeasure.is_invariant_under>` (R2 of
    #434; until 2026-09-03 an ``is_invariant`` method lived here and dragged
    the measure module into this one's imports).

    A frozen value type, like every tag it wraps and every
    :class:`~orpheus.numerics.manifold.Manifold` that carries one as a field:
    `[M]` 2026-09-02, with a writable slot, ``g._tag = …`` succeeded and moved
    ``hash(quotient)`` from under three memos keyed on it.
    """

    #: The tag — the group's IDENTITY and name: a named enum member or a
    #: parameterised family.  Every question about the group is answered by
    #: its :attr:`realization`, derived from the tag once and cached.
    _tag: "_NamedSubgroup | Cn | Dnh | Mirror | SO2 | O2"

    # Pre-instantiated named singletons — assigned at module scope below,
    # after the class and ``_NamedSubgroup`` are fully defined.  Declared
    # here as ClassVars so the public ``SubgroupOfO3.Trivial`` / ``.OctahedralOh``
    # / ... access surface is statically known. ``SO2``, ``O2`` and ``Mirror``
    # are NOT here: all three are axis-parameterised constructors below.
    Trivial: ClassVar["SubgroupOfO3"]
    Dinfh: ClassVar["SubgroupOfO3"]
    OctahedralOh: ClassVar["SubgroupOfO3"]
    IcosahedralIh: ClassVar["SubgroupOfO3"]
    SO3: ClassVar["SubgroupOfO3"]
    O3: ClassVar["SubgroupOfO3"]

    def __post_init__(self) -> None:
        # C_1 IS the trivial group: one group, one spelling — #432's naming
        # law applied to the finite families.  `[M]` 2026-09-02 a second
        # value for {e} gave a different orbit_stabiliser, a different hash,
        # and an EMPTY `_maximal` (each spelling "strictly" contained the
        # other).  Normalised here, on the type, so no constructor can mint
        # the second spelling.
        if isinstance(self._tag, Cn) and self._tag.n == 1:
            object.__setattr__(self, "_tag", _NamedSubgroup.Trivial)

    # --- Constructors ----------------------------------------------------

    @classmethod
    def _from_named(cls, name: _NamedSubgroup) -> "SubgroupOfO3":
        return cls(name)

    @classmethod
    def Cn(cls, n: int) -> "SubgroupOfO3":
        """Cyclic group :math:`C_n` of :math:`n`-fold proper rotations.
        ``Cn(1)`` IS :attr:`Trivial` — one group, one spelling."""
        return cls(Cn(n))

    @classmethod
    def Dnh(cls, n: int) -> "SubgroupOfO3":
        """Dihedral group :math:`D_{nh}` of order :math:`4n`."""
        return cls(Dnh(n))

    @classmethod
    def Mirror(cls, axis: str) -> "SubgroupOfO3":
        r"""Reflection group :math:`\{e, \sigma_a\}` about the plane
        normal to ``axis``. Replaced the parameter-free ``Z2`` on
        2026-08-02 — see :class:`Mirror` for why the plane is not a
        convention."""
        return cls(Mirror(axis))

    @classmethod
    def SO2(cls, axis: str) -> "SubgroupOfO3":
        r"""The axial rotation group :math:`C_\infty` about ``axis``.
        Replaced the parameter-free ``SO2`` entry on 2026-09-01 — see
        :class:`SO2` for why the axis is not a convention."""
        return cls(SO2(axis))

    @classmethod
    def O2(cls, axis: str) -> "SubgroupOfO3":
        r"""The stabiliser :math:`O(2)_a = C_{\infty v}` of ``axis`` — the
        rotations about it AND the mirrors through it. The group a
        :math:`\mu`-orbit space is named by (#432, 2026-09-02) — see
        :class:`O2`."""
        return cls(O2(axis))

    @property
    def rotation_axis(self) -> int | None:
        r"""The coordinate index of an AXIAL group's axis, else ``None``.

        ``SO2("x") → 0``, ``O2("x") → 0``, ``…("y") → 1``, ``…("z") → 2``:
        the two families whose every element fixes one coordinate axis
        pointwise — the proper rotations about it, and (for
        :math:`O(2)_a`, since 2026-09-02) those together with the mirrors
        through it. Every other subgroup answers ``None``, including the
        ones that CONTAIN axial rotations without fixing the axis
        (:math:`D_{\infty h}` flips it by :math:`\sigma_h`; :math:`SO(3)`
        and :math:`O(3)` move it) and the finite cyclic groups, which are
        realized about :math:`z` without naming it. The continuous dual of
        :attr:`mirror_axis`: it identifies the axis whose polar interval the
        group's orbit space on :math:`S^2` IS, which is what an orbit-space
        derivation (``manifold._sphere_mod_o2`` needs the invariant
        coordinate) and the embedding of a polar marginal
        (``invariance._embedded_nodes`` puts :math:`\mu` on THIS axis) read off it.
        """
        tag = self._tag
        if isinstance(tag, (SO2, O2)):
            return AXIS_INDEX[tag.axis]
        return None

    @property
    def is_trivial(self) -> bool:
        r"""``True`` for :math:`\{e\}` — the one group that acts trivially on
        every orbit space, admitted on every base by the catalogue door
        (which reads this, not the NAME — #434 R4)."""
        return self._tag is _NamedSubgroup.Trivial

    @property
    def realization(self) -> "Realization":
        r"""The group as what it IS — its identity component and one
        representative per connected component (:class:`Realization`),
        derived from the tag once (:func:`_realize`) and cached.  Every
        predicate on this class is a computation on it.  ⚠ Not to be confused
        with :attr:`~orpheus.numerics.manifold.Quotient.realization`, which
        is the MANIFOLD an orbit space is (``Ball(2)`` for
        :math:`S^2/\sigma_y`): a group's realization is its structure, an
        orbit space's is its shape."""
        return _realize(self._tag)

    @property
    def dim(self) -> int:
        r"""The group's dimension as a Lie group — 0 for a finite member, 1 for
        the axial families and :math:`D_{\infty h}`, 3 for :math:`SO(3)` and
        :math:`O(3)`.  Read off the identity component, which is the half that
        has one.  ⚠ **NOT** what an orbit space's dimension law subtracts: that
        is the generic ORBIT's dimension (:meth:`generic_orbit_dimension`), and
        ``dim H`` is wrong for :math:`O(3)` on :math:`S^2` (3, orbits 2) and
        :math:`SO(3)` on :math:`\mathbb{R}^3` — a law spelled ``base.dim -
        by.dim`` would give −1 and 0 (#434 R4, 2026-09-03)."""
        return self.realization.component.dim

    def generic_orbit_dimension(self, points: np.ndarray) -> int:
        r"""The dimension of this group's generic orbit through ``points`` —
        the largest :math:`G^0`-orbit over the probe set
        (:meth:`IdentityComponent.orbit_dimension`; the coset representatives
        move orbits without changing their dimension).  ⚠ NOT :attr:`dim`:
        :math:`O(3)` has :math:`\dim 3` and 2-dimensional orbits on
        :math:`S^2`.  The one surface :mod:`~orpheus.numerics.manifold`'s
        dimension law reads."""
        return self.realization.component.orbit_dimension(points)

    @property
    def orbit_stabiliser(self) -> "SubgroupOfO3":
        r"""The largest subgroup of :math:`O(3)` with this group's orbits —
        the group an orbit space is NAMED BY.

        Two groups with the same orbits on :math:`S^2` have the same
        invariant ring, hence ONE orbit space; so a catalogue entry
        :math:`S^2/H` is keyed on the LARGEST such group, and asking for
        the quotient by a smaller orbit-equivalent one is refused with this
        answer — at construction, by
        :class:`~orpheus.numerics.manifold.Quotient`, and at the catalogue
        door (#432, 2026-09-02).  Well-defined: the elements of
        :math:`O(3)` carrying every :math:`H`-orbit onto itself form a group
        containing :math:`H` with the same orbits, and every group with
        those orbits lies inside it.

        Exactly two FAMILIES are not their own stabiliser (`[M]` 4 of 26
        distinct members: :math:`SO(2)_x`, :math:`SO(2)_y`, :math:`SO(2)_z`,
        :math:`SO(3)`) — the CONNECTED groups, whose improper extension fixes
        every invariant:

        * :math:`SO(2)_a \to O(2)_a`:
          :math:`\mathbb{R}[x]^{SO(2)_a} = \mathbb{R}[x]^{O(2)_a} =
          \mathbb{R}[x_a,\, x_b^2 + x_c^2]` — a vertical mirror carries
          each constant-:math:`\mu` circle onto itself (:class:`O2`).
        * :math:`SO(3) \to O(3)`: :math:`\mathbb{R}[x]^{SO(3)} =
          \mathbb{R}[x]^{O(3)} = \mathbb{R}[\,|x|^2\,]` — the orbits are
          the spheres, and :math:`S^2/SO(3) = S^2/O(3)` is a point.

        Every other member answers itself.  `[R]` for a FINITE group an
        isometry that carries every orbit onto itself agrees with one
        element on an open set of points, hence everywhere — so
        :math:`\{e\}`, :math:`C_n`, :math:`D_{nh}` (with :math:`D_{1h}`),
        :math:`\sigma_a`, :math:`O_h`, :math:`I_h` are each their own;
        :math:`D_{\infty h}`, :math:`O(2)_a` and :math:`O(3)` are the
        stabilisers of their own orbit families (the latitude-circle pairs
        about :math:`z`, the circles about :math:`a`, the sphere).
        `[M]` 2026-09-02 the two non-identity cases are the two whose
        orbit-space entry the two groups would otherwise DUPLICATE:
        ``Quotient.descending_slots`` under ``O2(a)`` and under the
        pre-#432 ``SO2(a)`` entry are ``array_equal`` at every axis and
        :math:`L \le 4` — the invariant rings realized, not asserted.
        """
        realization = self.realization
        if realization.is_finite:
            return self
        if realization.component.dim == 3:
            return SubgroupOfO3.O3
        # A torus about a: the orbits are unions of circles about a.  When
        # the group acts trivially on the circles' heights (G ⊆ O(2)_a) the
        # partition's stabiliser is O(2)_a; when it pairs ±μ (D_∞h) it is the
        # stabiliser of the axis as a LINE, D_∞h(a) — the group itself.
        axial = SubgroupOfO3.O2(_letter_of(realization.component.axis))
        if axial.contains(self):
            return axial
        if self.contains(axial):
            return self
        raise NotImplementedError(
            f"the orbit stabiliser of {self.name} — a group with a torus "
            f"component neither inside nor containing {axial.name} — is not "
            f"spelled; no shipped member reaches this."
        )

    @property
    def mirror_axis(self) -> int | None:
        r"""The negated coordinate index for a single-mirror group, else None.

        ``Mirror("x") → 0``, ``Mirror("y") → 1``, ``Mirror("z") → 2``;
        every other subgroup (including those that CONTAIN mirrors)
        answers ``None`` — the property identifies the group whose sole
        non-identity element is one coordinate reflection, which is
        what a fold-aware consumer (the folded quadrature's σ-even
        harmonic sub-basis, Q5.6) needs to know.
        """
        tag = self._tag
        if isinstance(tag, Mirror):
            return AXIS_INDEX[tag.axis]
        return None

    # --- Repr / name -----------------------------------------------------

    def __repr__(self) -> str:
        tag = self._tag
        if isinstance(tag, _NamedSubgroup):
            return f"SubgroupOfO3.{tag.name}"
        if isinstance(tag, Cn):
            return f"SubgroupOfO3.Cn({tag.n})"
        if isinstance(tag, Dnh):
            return f"SubgroupOfO3.Dnh({tag.n})"
        if isinstance(tag, Mirror):
            return f"SubgroupOfO3.Mirror({tag.axis!r})"
        if isinstance(tag, SO2):
            return f"SubgroupOfO3.SO2({tag.axis!r})"
        if isinstance(tag, O2):
            return f"SubgroupOfO3.O2({tag.axis!r})"
        return f"SubgroupOfO3({tag!r})"

    @property
    def name(self) -> str:
        """Human-readable name for the group (used in error messages)."""
        tag = self._tag
        if isinstance(tag, _NamedSubgroup):
            return tag.value
        if isinstance(tag, Cn):
            return f"C_{tag.n}"
        if isinstance(tag, Dnh):
            return f"D_{tag.n}h"
        if isinstance(tag, Mirror):
            return f"sigma_{tag.axis}"
        if isinstance(tag, SO2):
            # The orbit-space catalogue keys on this string, so an
            # axis-blind name would merge three quotients into one entry.
            return f"SO2_{tag.axis}"
        if isinstance(tag, O2):
            # `S^2/O2_x` is the slab's orbit space's name (#432).
            return f"O2_{tag.axis}"
        return repr(tag)

    # --- Containment lattice --------------------------------------------

    def contains(self, other: "SubgroupOfO3") -> bool:
        r"""``True`` iff ``other ⊆ self`` in the :math:`O(3)` lattice.

        Implements :eq:`subgroup-of-o3-containment` by COMPUTATION on the two
        realizations (:meth:`Realization.contains`): the other group's
        identity component must lie in this one's — a torus in a torus iff the
        axes are parallel, anything in :math:`SO(3)`, :math:`\{e\}` in
        everything — and every one of the other group's coset representatives
        must be an ELEMENT of this group (:meth:`Realization.contains_element`,
        decided by which coset it lands in).  Containment is literal subgroup
        containment *in the standard setting* (principal axis along
        :math:`z`, vertex on :math:`x`), which is the question the
        quadrature-selection gate asks because a rule and its geometry are used
        in one frame; subconjugacy (:math:`\exists g: gHg^{-1} \le K`) is a
        different question and would get its own verb.

        The relations this computes, for the record and for the reader who
        wants the lattice at a glance: :math:`C_n \subseteq C_m` and
        :math:`D_{nh}` iff :math:`n \mid m`; :math:`C_n \subseteq SO(2)_z
        \subseteq D_{\infty h}` for every :math:`n`; :math:`D_{nh} \subseteq
        D_{\infty h}` and :math:`\not\subseteq O(2)_a` (its in-plane
        :math:`C_2` axes and :math:`\sigma_h` flip the axis); :math:`\sigma_b
        \subseteq O(2)_a` iff :math:`b \ne a`; :math:`D_{1h} = \{e, C_2^x,
        \sigma_y, \sigma_z\} \subseteq O(2)_x`; :math:`O_h \supseteq C_4,
        D_{4h}, C_2, D_{2h}` and the three coordinate mirrors (`[M]` — the
        comment that once called these "out of scope for the static lattice"
        was stale the day the finite arm was computed); :math:`O_h, I_h,
        D_{\infty h}, O(2)_a \not\subseteq SO(3)` (improper elements);
        :math:`SO(2)_a \subseteq O(2)_a \subseteq O(3)`; :math:`O(2)_z
        \subseteq D_{\infty h} = O(2)_z \times \{e, \sigma_h\}` and no
        other axis.  None of these is written down anywhere: each is what the
        two realizations say.

        Reflexivity (``self.contains(self) == True``) is a THEOREM of the
        realization — `[M]` 2026-09-03, ``realization.contains(itself)`` is
        ``True`` on all 31 expressible members — and the equality fast path
        in :func:`_tags_contain` is a cost memo on top of it (:math:`I_h`
        reflexivity is otherwise 120 × 120 element comparisons), not the
        reason it holds.

        Parameters
        ----------
        other : SubgroupOfO3
            The candidate subgroup.

        Returns
        -------
        bool
            ``True`` iff every element of ``other`` is also an element
            of ``self``.
        """
        return _tags_contain(self._tag, other._tag)

    # --- Invariance check -----------------------------------------------

    def is_normalised_by(
        self, motion: RigidMotion, *, atol: float = _ELEMENT_ATOL,
    ) -> bool:
        r"""``True`` iff ``motion`` normalises this group, :math:`g H g^{-1}
        = H` — the condition under which :math:`g` DESCENDS to the orbit
        space :math:`M/H`.

        An isometry carries :math:`H`-orbits onto :math:`H`-orbits exactly
        when it normalises :math:`H`; then :math:`[p] \mapsto [g\,p]` is a
        function of the orbit and :math:`g` acts on :math:`M/H`
        (:meth:`Quotient.induced_action
        <orpheus.numerics.manifold.Quotient.induced_action>`). Otherwise
        the image of an orbit depends on the representative and there is
        no action to speak of. Decided EXACTLY, family by family:

        * a FINITE group — the conjugated element set equals the element
          set (:meth:`RigidMotion.conjugated_by
          <orpheus.geometry.transformation.RigidMotion.conjugated_by>`);
        * :math:`SO(2)_a` and :math:`O(2)_a` — :math:`g` carries the line
          of :math:`a` onto itself, :math:`g\hat e_a = \pm\hat e_a`
          (the rotations about :math:`a` are conjugated to the rotations
          about :math:`g\hat e_a`, and the mirrors through :math:`a` to the
          mirrors through :math:`g\hat e_a`);
        * :math:`D_{\infty h}` — likewise for :math:`z`;
        * :math:`\{e\}`, :math:`SO(3)`, :math:`O(3)` — normal in
          :math:`O(3)`, so by everything.

        `[M]` 2026-09-02: :math:`C_4` about :math:`z` does not normalise
        :math:`\sigma_y` (it conjugates it to :math:`\sigma_x`), so it does
        not act on :math:`S^2/\sigma_y`; :math:`\sigma_x` does, and acts on
        the disk as :math:`(x, z) \mapsto (-x, z)`.

        Asked of the motion's LINEAR part: a point group acts on
        directions, and a translation does not move a direction — the
        convention :meth:`Quadrature.ordinate_permutation
        <orpheus.numerics.quadrature.directional.Quadrature.ordinate_permutation>`
        already states (a wrap's ordinate permutation is the identity), so
        a translated deck element normalises exactly what its rotation
        does.  Decided on the realization (:meth:`Realization.is_normalised_by`):
        :math:`g` must carry the Lie algebra onto itself (a torus: :math:`g\hat
        a = \pm\hat a`) and conjugate every coset representative into the
        group — which for a FINITE group is the conjugated element set equalling
        the element set, and for :math:`O(2)_a` / :math:`D_{\infty h}` the line
        criterion above, both from one body.
        """
        return self.realization.is_normalised_by(motion.linear_part, atol=atol)

    def normalises(
        self, other: "SubgroupOfO3", *, atol: float = _ELEMENT_ATOL,
    ) -> bool:
        r"""``True`` iff every element of ``self`` normalises ``other`` —
        ``self`` ACTS on the orbit space :math:`M/\text{other}`.

        Decided EXACTLY on the realization (:meth:`Realization.normalises`):
        the identity component through the Lie algebra — :math:`G^0
        \subseteq N(H)` iff :math:`[\mathfrak g, \mathfrak h] \subseteq
        \mathfrak h` and :math:`X - \mathrm{Ad}_s X \in \mathfrak h` for
        every generator :math:`X` and every coset representative :math:`s` of
        :math:`H` (:meth:`IdentityComponent.normalises` carries the proof) —
        and the coset representatives of ``self`` one by one.  On a finite
        ``other`` the Lie condition reads "every element commutes with the
        rotation generator :math:`[\hat e_a]_\times`" (conjugation by a
        connected group into a discrete one is constant); on :math:`O(2)_b` it
        reads :math:`a = b`; :math:`\{e\}`, :math:`SO(3)` and :math:`O(3)` are
        normal in everything.  Never sampled (ERR-072).
        """
        return self.realization.normalises(other.realization, atol=atol)

    @property
    def identity_component(self) -> "SubgroupOfO3":
        r"""The connected component of the identity, :math:`G^0` — as a NAMED
        member: :attr:`Trivial` for every finite member (a finite subgroup of
        :math:`O(3)` is discrete, so its identity component is :math:`\{e\}`),
        :math:`SO(2)_a` for the axial families and :math:`D_{\infty h}`,
        :math:`SO(3)` for :math:`SO(3)` and :math:`O(3)`.

        The half of a continuous group whose action on a finite point set is
        decided by CONNECTEDNESS rather than by enumeration: its orbits are
        connected, so it fixes every point of any finite invariant set
        (:meth:`IdentityComponent.fixes`).  `[M]` until 2026-09-02 this
        returned the group ITSELF for every finite member — wrong on 12 of 22
        expressible members, contradicting its own docstring on 11, and
        invisible because it had zero readers; the two places that needed the
        component destructured a private helper instead.  Now it is what the
        realization says, and the dimension law on
        :class:`~orpheus.numerics.manifold.Quotient` reads it.
        """
        component = self.realization.component
        if component.dim == 0:
            return SubgroupOfO3.Trivial
        if component.dim == 3:
            return SubgroupOfO3.SO3
        return SubgroupOfO3.SO2(_letter_of(component.axis))

    def generic_images(self, points: "np.ndarray") -> list["np.ndarray"]:
        r"""The images of ``points`` under a GENERIC set of this group's elements.

        The surface an orbit-space entry's isotypic probe needs
        (:meth:`Quotient.descending_slots
        <orpheus.numerics.manifold.Quotient.descending_slots>`, #429 tracker
        3.4): a function is constant on the group's orbits iff it agrees at
        every point with its value at every image. For a FINITE group the
        generic set is every element (the memoised closure of the
        realization). For :math:`SO(2)_a` it is rotations about :math:`a` by
        INCOMMENSURATE angles, and for :math:`O(2)_a` those rotations
        together with each composed with one vertical mirror — a generic
        element of each component: a finite sample of a continuous group
        generates a finite SUBgroup, and a sample of right angles generates
        :math:`C_4`, which `[M]` 2026-09-02 falsely admits the
        :math:`m = \pm 4` real-harmonic slots at :math:`L \ge 4`
        (``vv-principles`` #13 — this is the trap
        :func:`~orpheus.numerics.invariance.is_invariant_under` avoids
        by deciding continuous groups EXACTLY; a probe of FUNCTIONS cannot,
        so it samples where no finite subgroup can hide). :math:`D_{\infty h}`
        has an axis (:math:`z`) and four components, so since 2026-09-02 (R1
        of #434) its realization samples the torus composed with each coset
        representative — `[M]` 6 angles × 4 components = 24 images, where the
        retired tag dispatch refused it as "continuous". :math:`SO(3)` and
        :math:`O(3)` have no axis to rotate about and still refuse until a
        consumer needs them (:meth:`Realization.generic_images`).

        ``points`` is ``(N, 3)``; each image is ``(N, 3)``, the linear part
        of the element applied row-wise (every realized element fixes the
        origin).
        """
        return self.realization.generic_images(points)


# ---------------------------------------------------------------------------
# The realization — a closed subgroup of O(3) as what it IS
# ---------------------------------------------------------------------------


def _realized_ops(tag) -> "list[RigidMotion] | None":
    r"""The FINITE generating set of ``tag``, or ``None`` if continuous.

    Every realization is in the **standard setting**: principal axis along
    z, vertex on the x-axis.  These are GENERATING sets, not groups
    (:math:`D_{nh}` ships :math:`2n+1` matrices for a group of order
    :math:`4n`); :func:`_realize` closes them.  Kept as the finite half's
    entry point because the exactness gates read it by name.
    """
    if isinstance(tag, _NamedSubgroup):
        if tag is _NamedSubgroup.Trivial:
            return [RigidMotion.identity(3)]
        if tag is _NamedSubgroup.OctahedralOh:
            return _octahedral_ops()
        if tag is _NamedSubgroup.IcosahedralIh:
            return _icosahedral_ops()
        return None  # Dinfh / SO3 / O3 are continuous — no finite set
    if isinstance(tag, (SO2, O2)):
        return None  # continuous, whatever the axis
    if isinstance(tag, Mirror):
        # One generator, so the closed group is {e, sigma_a}: order 2.
        return _reflections(tag.axis)
    if isinstance(tag, Cn):
        return _cyclic_ops(tag.n)
    if isinstance(tag, Dnh):
        return _cyclic_ops(tag.n) + _reflections("z") + _vertical_mirrors(tag.n)
    raise ValueError(
        f"{tag!r} is not a member of the lattice; a tag with no arm here "
        f"would otherwise read as CONTINUOUS, which is the bare-return "
        f"defect the retired `_contains` shipped."
    )


def _close_group(ops: "list[RigidMotion]") -> "list[RigidMotion]":
    """The group generated by ``ops`` — closure under composition.

    A thin wrapper on :func:`~orpheus.geometry.transformation.close_group`,
    kept because callers here always mean *the point group in the standard
    3-D setting*. The rounding-key trick, the quadratic-not-cubic membership
    test and the finite-group cap all live in the core, where the cap
    raises a NAMED refusal rather than a bare ``ValueError``.
    """
    return list(_close_rigid_motions(ops, dimension=3))


def _skew(vector: np.ndarray) -> np.ndarray:
    r"""The infinitesimal rotation :math:`[v]_\times` about ``vector`` — the
    generator with :math:`R_v(\theta) = \exp(\theta\,[\hat v]_\times)`."""
    v = np.asarray(vector, dtype=float)
    return np.array(
        [[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]]
    )


def _axis_of(generator: np.ndarray) -> np.ndarray:
    """The unit axis of a one-dimensional generator — :func:`_skew`'s inverse."""
    axis = np.array([generator[2, 1], generator[0, 2], generator[1, 0]])
    return axis / np.linalg.norm(axis)


def _letter_of(axis: np.ndarray) -> str:
    """The coordinate letter of a unit axis along a coordinate direction.

    The realization carries AXES as vectors; the named members carry them as
    letters.  A group about a non-coordinate axis has no named member yet,
    and this is where that refusal lives.
    """
    for i, letter in enumerate(AXIS_LETTER):
        if abs(abs(float(axis[i])) - 1.0) <= _ELEMENT_ATOL:
            return letter
    raise NotImplementedError(
        f"no named member of the lattice is about the axis {axis}; the named "
        f"axial families are about x, y and z."
    )


def _perpendicular_letter(axis: str) -> str:
    r"""The first coordinate axis other than ``axis`` — the normal of the
    vertical mirror representing :math:`O(2)_a`'s improper coset.  Any plane
    containing the axis generates that coset; a coordinate mirror is a
    bit-exact signed diagonal."""
    return next(b for b in AXIS_LETTER if b != axis)


def _in_span(generators: "tuple[np.ndarray, ...]", matrix: np.ndarray, atol: float) -> bool:
    r"""``matrix ∈ span(generators)`` inside :math:`\mathfrak{so}(3)`, to ``atol``.

    Three cases by dimension — the only three a Lie subalgebra of
    :math:`\mathfrak{so}(3)` has: the zero algebra (the matrix must vanish),
    a line (the component orthogonal to the generator must vanish), the
    whole algebra (every skew matrix belongs).
    """
    if len(generators) == 0:
        return bool(np.abs(matrix).max() <= atol)
    if len(generators) == 3:
        return True
    (generator,) = generators
    coefficient = float(np.sum(matrix * generator) / np.sum(generator * generator))
    return bool(np.abs(matrix - coefficient * generator).max() <= atol)


@dataclass(frozen=True, eq=False)
class IdentityComponent:
    r"""The connected component of the identity, :math:`G^0 = \exp\mathfrak g`.

    ``generators`` is a basis of :math:`\mathfrak g \subseteq
    \mathfrak{so}(3)`: empty for a finite (discrete) group, one
    :math:`[\hat a]_\times` for a torus about :math:`\hat a`, three for
    :math:`SO(3)` — and nothing else, because :math:`\mathfrak{so}(3)` has no
    two-dimensional subalgebra.  Every exact statement about a continuous
    group in this module is a statement about this half plus finitely many
    coset representatives (:class:`Realization`), never about a sample of
    its elements (ERR-072): a connected group's action on a FINITE set is
    decided by its Lie algebra, because its orbits are connected.

    Not compared by value (``eq=False``): a component is derived from a tag,
    and the tag is the identity.
    """

    generators: "tuple[np.ndarray, ...]"

    def __post_init__(self) -> None:
        # `dim` is a RANK, and every predicate below spells it `len(generators)`.
        # so(3) has subalgebras of dimension 0, 1 and 3 only, so a dependent or
        # non-skew basis does not merely mis-answer: at length 3 `_in_span`
        # returns True unconditionally and a LINE claims to be so(3).  `[M]`
        # 2026-09-03 (elegance review) all three states were constructible;
        # 0 of the 31 shipped members trip the guard.
        n = len(self.generators)
        if n not in (0, 1, 3):
            raise ValueError(
                f"a Lie subalgebra of so(3) has dimension 0, 1 or 3, not {n}."
            )
        for X in self.generators:
            if X.shape != (3, 3) or np.abs(X + X.T).max() > _ELEMENT_ATOL:
                raise ValueError(
                    f"a generator must be a 3x3 skew-symmetric matrix, got\n{X}"
                )
        if n and np.linalg.matrix_rank(
            np.stack([X.ravel() for X in self.generators]), tol=_ELEMENT_ATOL,
        ) != n:
            raise ValueError(
                f"the {n} generators of this component are linearly dependent, "
                f"so its dimension is not {n}; `dim` is a rank."
            )

    @property
    def dim(self) -> int:
        return len(self.generators)

    @property
    def axis(self) -> np.ndarray:
        """The torus axis — a unit vector.  Only a one-dimensional component
        has one; asking :math:`\\{e\\}` or :math:`SO(3)` is refused."""
        if self.dim != 1:
            raise ValueError(
                f"only a torus has an axis; this component has dimension {self.dim}."
            )
        return _axis_of(self.generators[0])

    def contains_element(self, motion: RigidMotion, *, atol: float = _ELEMENT_ATOL) -> bool:
        r"""``motion`` :math:`\in G^0`.  By dimension: :math:`\{e\}` admits the
        identity alone; a torus admits a PROPER motion fixing its axis;
        :math:`SO(3)` admits every proper motion."""
        return self._contains_matrix(motion.linear, atol)

    def _contains_matrix(self, Q: np.ndarray, atol: float) -> bool:
        r"""The membership test on the orthogonal MATRIX — the form a coset
        test composes (:meth:`Realization.contains_element`), so that asking
        whether :math:`r^{-1} g` lies in :math:`G^0` costs one matrix product
        and not the minting of an element (with its orthogonality check)."""
        if self.dim == 0:
            return bool(np.abs(Q - _IDENTITY_3).max() <= atol)
        if np.linalg.det(Q) < 0.0:
            return False
        if self.dim == 3:
            return True
        axis = self.axis
        return bool(np.abs(Q @ axis - axis).max() <= atol)

    def contains(self, other: "IdentityComponent", *, atol: float = _ELEMENT_ATOL) -> bool:
        r""":math:`\mathfrak h \subseteq \mathfrak g` — a torus lies in a torus
        iff the axes are parallel, anything lies in :math:`SO(3)`,
        :math:`\{e\}` lies in everything."""
        return all(_in_span(self.generators, Y, atol) for Y in other.generators)

    def is_normalised_by(self, motion: RigidMotion, *, atol: float = _ELEMENT_ATOL) -> bool:
        r""":math:`\mathrm{Ad}_g\,\mathfrak g = \mathfrak g` — for a torus,
        :math:`g` carries the axis LINE onto itself, :math:`g\hat a = \pm\hat
        a` (:math:`Q[\hat a]_\times Q^{\mathsf T} = \det Q\,[Q\hat a]_\times`,
        and the sign is invisible to a span)."""
        Q = motion.linear
        return all(_in_span(self.generators, Q @ X @ Q.T, atol) for X in self.generators)

    def normalises(self, other: "Realization", *, atol: float = _ELEMENT_ATOL) -> bool:
        r"""Does :math:`G^0` normalise :math:`H` — EXACTLY, through the Lie algebra.

        :math:`G^0 \subseteq N(H)` iff :math:`\mathfrak g \subseteq
        \mathrm{Lie}\,N(H)`, and

        .. math::

           \mathrm{Lie}\,N(H) \;=\; \{\,X \in \mathfrak{so}(3) :
           [X, \mathfrak h] \subseteq \mathfrak h \ \text{ and }\
           X - \mathrm{Ad}_s X \in \mathfrak h \ \ \forall\, s \in H \,\}.

        Necessity is the derivative at :math:`t = 0` of :math:`\exp(tX)\,s\,
        \exp(-tX)\,s^{-1} \in H^0`.  Sufficiency: with :math:`Y_s = X -
        \mathrm{Ad}_s X \in \mathfrak h`, the curve :math:`f(t) = \exp(tX)\,s\,
        \exp(-tX)\,s^{-1}` satisfies :math:`f'(t) f(t)^{-1} =
        \mathrm{Ad}_{\exp(tX)} Y_s \in \mathrm{Ad}_{\exp(tX)}\mathfrak h =
        \mathfrak h` (the first condition makes :math:`\mathfrak h`
        :math:`\mathrm{Ad}_{\exp(tX)}`-stable), so :math:`f` stays in
        :math:`H^0` — and conjugation by :math:`\exp(tX)` carries each
        component of :math:`H` onto itself.  Checking the second condition
        on the coset REPRESENTATIVES alone suffices for all of :math:`H`:
        for :math:`s\,h` with :math:`h = \exp Y`, :math:`Y \in \mathfrak h`,
        :math:`X - \mathrm{Ad}_{sh}X = (X - \mathrm{Ad}_s X) + \mathrm{Ad}_s
        (X - \mathrm{Ad}_h X)`, where :math:`X - \mathrm{Ad}_h X = -(\mathrm{ad}_Y
        X + \tfrac12 \mathrm{ad}_Y^2 X + \dots) \in \mathfrak h` by the first
        condition, and :math:`\mathrm{Ad}_s` preserves :math:`\mathfrak h`
        because :math:`H^0` is normal in :math:`H`.

        On a FINITE :math:`H` (:math:`\mathfrak h = 0`) this reads "every
        element commutes with :math:`X`": conjugation by a connected group
        into a discrete one is constant.  On :math:`O(2)_b` it reads
        :math:`\hat a \parallel \hat b`, since :math:`[[\hat a]_\times,
        [\hat b]_\times] = [\hat a \times \hat b]_\times` lies in the line of
        :math:`[\hat b]_\times` only when the cross product vanishes, and
        then :math:`X - \mathrm{Ad}_{\sigma_v} X = 2X \in \mathfrak h` (an
        improper :math:`g` sends :math:`[v]_\times` to :math:`-[gv]_\times`).
        On :math:`D_{\infty h}` it reads :math:`\hat a = \hat z`.  One body,
        every answer the retired per-family arms gave.
        """
        h = other.component
        for X in self.generators:
            if not all(_in_span(h.generators, X @ Y - Y @ X, atol) for Y in h.generators):
                return False
            for s in other.representatives:
                Q = s.linear
                if not _in_span(h.generators, X - Q @ X @ Q.T, atol):
                    return False
        return True

    def fixes(self, points: np.ndarray, *, atol: float) -> bool:
        r"""``True`` iff :math:`G^0` fixes every one of ``points`` —
        :math:`Xp = 0` for every generator :math:`X` and every point.

        The exact criterion for a CONNECTED group's invariance of a finite
        set: its orbits are connected, and a connected orbit inside a finite
        set is a point.  For a torus about :math:`\hat a`, :math:`[\hat
        a]_\times p = \hat a \times p` vanishes iff :math:`p` lies ON the axis
        (:math:`|\hat a \times p|` is its distance from the axis); for
        :math:`SO(3)` iff :math:`p = 0`.

        This is the criterion that replaced the sampled :math:`\{0, 90, 180,
        270\}^\circ` check, which tested closure under :math:`C_4` and
        called it :math:`SO(2)` (ERR-072) — unsound in the dangerous
        direction, certifying non-invariant rules as a function of ``n_phi
        mod 4``.  Consequence, and it is intended: **no real angular
        cubature is SO(2)-invariant.**  A consumer needing "this rule
        respects a continuous azimuthal symmetry" is asking about the rule's
        exactness space, not about node-set closure.  The measures that DO
        pass are the polar marginals — a rule on :math:`S^2/O(2)_a` embeds
        ON axis :math:`a`, where the group acts trivially — and they pass for
        their own axis only.  The icosahedral-sample check that once stood
        in for :math:`SO(3)` had the same defect, twice over: :math:`-I \in
        I_h`, so the ``SO3`` and ``O3`` branches ran one operator set, 60 of
        whose 120 matrices are not in :math:`SO(3)` at all.
        """
        pts = np.asarray(points, dtype=float)
        if pts.ndim == 1:
            pts = pts[None, :]
        if pts.size == 0:
            return True
        return all(bool(np.abs(pts @ X.T).max() <= atol) for X in self.generators)

    def orbit_dimension(self, points: np.ndarray, *, atol: float = _ELEMENT_ATOL) -> int:
        r"""The dimension of the LARGEST :math:`G^0`-orbit through ``points`` —
        :math:`\max_p \operatorname{rank}[X p : X \in \mathfrak g]`, the
        tangent space of the orbit through :math:`p`; orbit dimension is upper
        semicontinuous, so the maximum over a generic probe set IS the generic
        value.  :meth:`fixes` is this at zero: ``fixes(pts)`` iff
        ``orbit_dimension(pts) == 0``.  The rank is taken at the element band,
        the same tolerance that validates this component's own ``dim``
        (:meth:`__post_init__`), so the two rank questions cannot answer
        differently.  What an orbit space's dimension law subtracts
        (:class:`~orpheus.numerics.manifold.Quotient`): `[D]` a torus on
        :math:`S^2` has 1-dimensional orbits, :math:`SO(3)` 2-dimensional ones
        (the sphere is one orbit), :math:`SO(3)` on :math:`\mathbb{R}^3` also
        2 (so :math:`\mathbb{R}^3/SO(3)` is one-dimensional).
        """
        if not self.generators:
            return 0
        pts = np.atleast_2d(np.asarray(points, dtype=float))
        return max(
            int(np.linalg.matrix_rank(np.stack([X @ p for X in self.generators]), tol=atol))
            for p in pts
        )



@dataclass(frozen=True, eq=False)
class Realization:
    r"""A closed subgroup of :math:`O(3)` as :math:`G = \bigsqcup_r r\,G^0` —
    its identity component and ONE representative per connected component,
    the identity first.  For a finite group the representatives are every
    element.

    Every predicate is one body: an element lies in :math:`G` iff it lies in
    some coset (:meth:`contains_element`); :math:`H \subseteq G` iff
    :math:`\mathfrak h \subseteq \mathfrak g` and every representative of
    :math:`H` is an element (:meth:`contains`); :math:`g` normalises
    :math:`G` iff it preserves :math:`\mathfrak g` and conjugates every
    representative into :math:`G` (:meth:`is_normalised_by`); :math:`G`
    normalises :math:`H` iff its component does, exactly, and each of its
    representatives does (:meth:`normalises`).  A hand-maintained relation
    table is a claim with no construction behind it, and this module shipped
    two such claims that were false; a realization cannot disagree with
    itself.
    """

    component: IdentityComponent
    representatives: "tuple[RigidMotion, ...]"

    def __post_init__(self) -> None:
        # G = ⊔ r·G⁰ has at least the identity coset, and "identity first" is
        # a contract `elements` and `generic_images` hand out in order.  `[M]`
        # 2026-09-03 (elegance review) the empty tuple constructed a "group"
        # whose contains_element(e) was False.
        if not self.representatives:
            raise ValueError(
                "a subgroup has at least the identity coset; `representatives` "
                "must not be empty."
            )
        if not self.representatives[0].approximately_equals(
            RigidMotion.identity(3), atol=_ELEMENT_ATOL,
        ):
            raise ValueError(
                "the first coset representative must be the identity — "
                "`elements` and `generic_images` hand the tuple out in that order."
            )

    @property
    def is_finite(self) -> bool:
        return self.component.dim == 0

    @property
    def elements(self) -> "tuple[RigidMotion, ...]":
        """Every element — defined for a FINITE group, where the
        representatives are the elements."""
        if not self.is_finite:
            raise ValueError(
                "a continuous group has no finite element set; ask its "
                "identity component and coset representatives instead."
            )
        return self.representatives

    def contains_element(self, motion: RigidMotion, *, atol: float = _ELEMENT_ATOL) -> bool:
        r"""``motion`` :math:`\in G`: it lies in some coset :math:`r\,G^0`,
        i.e. :math:`r^{-1} g \in G^0` for some representative :math:`r`."""
        Q = motion.linear
        return any(
            self.component._contains_matrix(r.linear.T @ Q, atol)
            for r in self.representatives
        )

    def contains(self, other: "Realization", *, atol: float = _ELEMENT_ATOL) -> bool:
        r""":math:`H \subseteq G`: :math:`H^0 \subseteq G^0` (connected, so
        through the Lie algebras) and every representative of :math:`H` is
        an element of :math:`G`."""
        return self.component.contains(other.component, atol=atol) and all(
            self.contains_element(s, atol=atol) for s in other.representatives
        )

    def is_normalised_by(self, motion: RigidMotion, *, atol: float = _ELEMENT_ATOL) -> bool:
        r""":math:`g G g^{-1} = G`: :math:`g` preserves :math:`\mathfrak g` and
        conjugates every representative into :math:`G` — for a finite group,
        the conjugated element set equals the element set."""
        return self.component.is_normalised_by(motion, atol=atol) and all(
            self.contains_element(r.conjugated_by(motion), atol=atol)
            for r in self.representatives
        )

    def normalises(self, other: "Realization", *, atol: float = _ELEMENT_ATOL) -> bool:
        r"""Every element of :math:`G` normalises :math:`H`: the identity
        component exactly (:meth:`IdentityComponent.normalises`), the
        representatives one by one."""
        return self.component.normalises(other, atol=atol) and all(
            other.is_normalised_by(r, atol=atol) for r in self.representatives
        )

    def generic_images(self, points: np.ndarray) -> "list[np.ndarray]":
        r"""The images of ``points`` under a GENERIC set of the group's elements.

        A finite group: every element.  A torus component: rotations by
        INCOMMENSURATE angles (:data:`_INCOMMENSURATE_ANGLES`), each composed
        with every coset representative — a generic element of each
        component, because a finite sample of a continuous group generates
        a finite SUBgroup, and a sample of right angles generates
        :math:`C_4`, which `[M]` 2026-09-02 falsely admits the :math:`m =
        \pm 4` real-harmonic slots at :math:`L \ge 4` (``vv-principles``
        #13; a probe of FUNCTIONS cannot be decided exactly the way a point
        set can, so it samples where no finite subgroup can hide).
        :math:`SO(3)` and :math:`O(3)` refuse until a consumer needs them.
        """
        pts = np.asarray(points, dtype=float)
        if self.component.dim == 3:
            raise NotImplementedError(
                "a group with SO(3) as its identity component has no rotation "
                "axis to sample about, and its generic images are not sampled "
                "another way; no consumer has needed them yet."
            )
        if self.is_finite:
            return [pts @ r.linear.T for r in self.representatives]
        rotations = [
            RigidMotion.rotation_about_axis(axis=self.component.axis, angle=theta)
            for theta in _INCOMMENSURATE_ANGLES
        ]
        return [
            pts @ (rotation @ r).linear.T
            for r in self.representatives
            for rotation in rotations
        ]


@functools.cache
def _tags_contain(outer, inner) -> bool:
    """:meth:`SubgroupOfO3.contains` on the two TAGS — memoised, because a
    containment answer is a pure function of them and the lattice walk asks
    :math:`O(n^3)` questions per measure, most of them repeats (`[M]` 420
    per walk on a slab rule — 203 distinct pairs, 217 cache hits,
    ``cache_info`` after one ``gauss_legendre(8).measure.symmetry_groups()``,
    2026-09-03)."""
    if outer == inner:
        return True
    return _realize(outer).contains(_realize(inner))


@functools.cache
def _realize(tag) -> Realization:
    r"""The realization of a tag — the ONE place the tag is read for structure.

    Continuous members by their Lie algebra and coset representatives:
    :math:`SO(2)_a = \exp(\mathbb R[\hat a]_\times)`; :math:`O(2)_a =
    SO(2)_a \sqcup SO(2)_a\,\sigma_v`; :math:`D_{\infty h} = SO(2)_z \sqcup
    SO(2)_z\sigma_h \sqcup SO(2)_z\sigma_v \sqcup SO(2)_z\sigma_h\sigma_v`;
    :math:`SO(3) = \exp\mathfrak{so}(3)`; :math:`O(3) = SO(3) \sqcup SO(3)(-I)`.
    Finite members by closing their generating set (:func:`_realized_ops`).
    Memoised on the tag — an immutable value — because the lattice walk asks
    :math:`O(n^3)` questions per measure and every one would otherwise re-close
    the same group (`[M]` a single walk once rebuilt :math:`I_h` 41 times,
    9.3 s of a 9.4 s walk).
    """
    identity = RigidMotion.identity(3)
    if isinstance(tag, (SO2, O2)):
        component = IdentityComponent((_skew(_axis_vector(tag.axis)),))
        if isinstance(tag, SO2):
            return Realization(component, (identity,))
        sigma_v = RigidMotion.reflection(
            normal=_axis_vector(_perpendicular_letter(tag.axis))
        )
        return Realization(component, (identity, sigma_v))
    if tag is _NamedSubgroup.Dinfh:
        sigma_h = RigidMotion.reflection(normal=_axis_vector("z"))
        sigma_v = RigidMotion.reflection(normal=_axis_vector("x"))
        return Realization(
            IdentityComponent((_skew(_axis_vector("z")),)),
            (identity, sigma_h, sigma_v, sigma_h @ sigma_v),
        )
    if tag is _NamedSubgroup.SO3 or tag is _NamedSubgroup.O3:
        so3 = IdentityComponent(tuple(_skew(_axis_vector(a)) for a in AXIS_LETTER))
        representatives = (
            (identity,) if tag is _NamedSubgroup.SO3
            else (identity, RigidMotion.inversion(3))
        )
        return Realization(so3, representatives)
    generators = _realized_ops(tag)
    if generators is None:  # a continuous tag every arm above should have taken
        raise ValueError(f"{tag!r} is continuous but has no realization arm")
    return Realization(IdentityComponent(()), tuple(_close_group(generators)))


def _group_elements(tag) -> "list[RigidMotion] | None":
    """Every element of a FINITE tag's group (memoised through
    :func:`_realize`), or ``None`` for a continuous one."""
    realization = _realize(tag)
    return list(realization.elements) if realization.is_finite else None


_INCOMMENSURATE_ANGLES: tuple[float, ...] = (
    1.0, float(np.sqrt(2.0)), float(np.e), 2.5, float(np.sqrt(7.0)),
    float(np.pi / 3.0 + 0.1),
)


def _axis_vector(axis: str) -> np.ndarray:
    """The unit vector along a named coordinate axis."""
    try:
        return np.eye(3)[AXIS_INDEX[axis]]
    except KeyError:
        raise ValueError(f"axis must be x/y/z, got {axis!r}") from None


def _reflections(axis: str) -> list[RigidMotion]:
    r"""The single coordinate mirror whose NORMAL is ``axis``.

    :math:`\sigma_x` is the reflection in the plane :math:`x = 0`, i.e. the
    one with normal :math:`\hat e_x` — not "a reflection through the x axis",
    which in :math:`\mathbb{R}^3` would be :math:`\mathrm{diag}(1,-1,-1)`, a
    half-turn with :math:`\det = +1`. The core makes the wrong reading
    unspellable by taking ``normal`` keyword-only; this wrapper exists to
    translate the tag's axis LETTER into that normal.
    """
    return [RigidMotion.reflection(normal=_axis_vector(axis))]


def _inversion_op() -> RigidMotion:
    r"""Inversion through the origin, :math:`x \to -x`.

    :math:`\det = (-1)^3 = -1`, and it fixes only the origin — so in
    :math:`\mathbb{R}^3` it is neither a reflection (which would fix a plane)
    nor a rotation (which would fix a line). That is why the core builds it
    from its definition rather than from either family.
    """
    return RigidMotion.inversion(3)


def _cyclic_ops(n: int) -> list[RigidMotion]:
    r""":math:`C_n`: the ``n`` proper rotations about the z-axis.

    The angles are the :math:`n`-th roots of unity, taken from
    :func:`~orpheus.numerics.roots_of_unity.roots_of_unity` as exact points
    on the circle rather than recomputed from :math:`\theta_k = 2\pi k/n`
    through :func:`numpy.cos`. Two reasons, and the second is the
    structural one:

    1. The symmetric angles have exactly-representable components
       (:math:`\cos(\pi/2) = 0`) that a rounded :math:`\pi/2` does not
       reproduce. `[M]` at :math:`n = 4` this takes the operator matrices
       from 30/36 exactly-representable entries to **36/36**, and the
       landing residual on ``product(4, 4)`` from ``2.1e-16`` to exactly
       ``0.0``.
    2. It is the SAME generator the azimuthal rules are built from
       (:func:`~orpheus.numerics.quadrature.periodic_trapezoid`). A checker
       that derives its own cosines is a second spelling of the very nodes
       it is checking — so the two sides could drift, and the tolerance
       absorbing that drift would be reported as the rule's own accuracy.
    """
    cos, sin = roots_of_unity(np.arange(n), n)
    return [
        RigidMotion.rotation_from_circle_point(
            plane=(_axis_vector("x"), _axis_vector("y")), point=(cos[k], sin[k])
        )
        for k in range(n)
    ]


def _vertical_mirrors(n: int) -> list[RigidMotion]:
    r"""The n vertical mirror planes of :math:`D_{nh}` (containing the z-axis).

    Built as the **coset** :math:`C_n\,\sigma_0` — that is,
    :math:`\sigma_k = R(2\pi k/n)\,\sigma_0`, where :math:`\sigma_0` is the
    :math:`\varphi = 0` mirror (normal :math:`\hat e_y`). This is the group
    fact :math:`D_n = C_n \sqcup C_n\sigma` spelled directly, and it has
    three consequences worth naming:

    - The mirror set IS the rotation set carried by one reflection, not a
      second trigonometric construction of the same angles. There is one
      source of truth for the :math:`n`-fold angles in this module.
    - :math:`\sigma_0` is a coordinate mirror, hence a bit-exact signed
      diagonal, so composing with it is a column sign flip and introduces
      **no** round-off. `[M]` the mirrors' landing residual is therefore
      exactly EQUAL to the rotations' on every rule measured — the mirror
      half can no longer be the worse half, which it previously was by
      1.7--4.8x because it went through a half-angle and a normalisation.
    - The half-angle disappears. The normal at :math:`k\pi/n + \pi/2` is
      a root of unity of order :math:`4n`, i.e. a *different* generator
      from the rule's own order-:math:`n` one; the coset form needs only
      the order-:math:`n` roots the rules themselves use.

    The plane placement is the convention that matters, and it is fixed
    by the standard :math:`D_{nh}` setting: the principal axis along z
    with a vertex on the **x-axis**, so :math:`\sigma_v` through
    :math:`\varphi = 0` is a symmetry. This is the setting every
    azimuthal rule in the tree is built in —
    :func:`~orpheus.numerics.quadrature.periodic_trapezoid` puts a node
    at :math:`\varphi = 0` for the unshifted (``NODE_ALIGNED``) family.
    Composing on the RIGHT by :math:`\sigma_0` is what places the k-th
    plane at :math:`k\pi/n`: :math:`\sigma_k` sends
    :math:`\varphi \mapsto 2\pi k/n - \varphi`, whose fixed azimuths are
    :math:`k\pi/n` and :math:`k\pi/n + \pi`.

    Previously the k-th **normal** was placed at :math:`k\pi/n`, i.e.
    every plane was rotated by :math:`\pi/2`. For even ``n`` that maps
    the normal set onto itself (:math:`\pi/2` is then a multiple of
    :math:`\pi/n`) and is invisible; for **odd** ``n`` it is a
    genuinely different set of planes, and
    ``Dnh(n).is_invariant(product(4, n))`` read ``False`` for
    ``n = 3, 5, 7`` while the rule is demonstrably closed under the
    :math:`\varphi = 0` mirror at every ``n``. Orthogonality,
    determinant, closure and group order are all preserved by a rotated
    mirror set, so none of those checks can see it. Two gates can, by
    different mechanisms: the plane-placement gate in
    ``tests/numerics/test_symmetry.py`` compares the placement against
    the setting, and the mirror-accuracy gate in
    ``tests/numerics/test_symmetry_exactness.py`` sees it as a RESIDUAL
    — rotating the planes takes the mirrors off the odd-``n`` rule
    entirely, so the landing blows up from ``2.5e-16`` to ``1.97e-01``.
    `[M]` that second gate reds on odd ``n`` and stays green on even
    ``n``, which is this paragraph's own claim, measured.
    """
    sigma_0 = RigidMotion.reflection(normal=_axis_vector("y"))
    return [rotation @ sigma_0 for rotation in _cyclic_ops(n)]


@functools.cache
def _octahedral_ops() -> list[RigidMotion]:
    """Generator set for the full octahedral group :math:`O_h` (48 elements).

    Memoised (2026-09-01): it is a constant, and until the realization
    was cached (:func:`_realize`, 2026-09-02) every containment question
    involving :math:`O_h` asked ``_realized_ops`` for it — see
    :func:`_icosahedral_ops` for the measurement.

    The full group is generated by:

    - All sign flips on the three coordinate axes (8 elements,
      :math:`(\\pm 1)^3`).
    - All coordinate permutations of :math:`(x, y, z)` (6 elements,
      :math:`S_3`).

    The Cartesian product gives :math:`8 \\times 6 = 48` matrices —
    exactly :math:`|O_h|`. The orbit-closure check applies all 48; a
    measure invariant under all of them is :math:`O_h`-invariant.
    """
    from itertools import product, permutations
    return [
        RigidMotion.signed_permutation(permutation=perm, signs=signs)
        for signs in product([-1.0, 1.0], repeat=3)
        for perm in permutations([0, 1, 2])
    ]


@functools.cache
def _icosahedral_ops() -> list[RigidMotion]:
    """Generator set for the icosahedral group :math:`I_h` (120 elements).

    ⚠ Memoised (2026-09-01), and it matters: this "generator set" is
    itself a full 120-element CLOSURE (~0.2 s).  `[M]` before the
    realization was cached (:func:`_realize`, 2026-09-02) a single
    ``symmetry_groups`` walk (then ``maximal_invariance_groups``, in this
    module) over a product rule rebuilt it
    **41 times** (9.3 s of a 9.4 s walk) once the axial family offered
    three axes instead of one — the retired containment dispatch asked
    for the generating set on every question, before any element cache
    was consulted.

    The icosahedral group is generated by a 5-fold rotation about a
    vertex axis and a 3-fold rotation about a face axis, plus
    inversion. We construct a representative orbit by:

    1. Listing the 12 vertices of a regular icosahedron.
    2. Generating proper rotations that map the first vertex to each
       other vertex (12 cosets) — gives the 60 elements of :math:`I`.
    3. Multiplying by inversion to extend to :math:`I_h` (120 elements).

    This is the standard construction; see Hamermesh 1962 §9.4.
    """
    # Icosahedron vertices (golden-ratio convention, normalised).
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    raw = np.array([
        [0,  1,  phi], [0,  1, -phi], [0, -1,  phi], [0, -1, -phi],
        [1,  phi, 0], [1, -phi, 0], [-1,  phi, 0], [-1, -phi, 0],
        [phi, 0,  1], [phi, 0, -1], [-phi, 0,  1], [-phi, 0, -1],
    ], dtype=float)
    raw /= np.linalg.norm(raw[0])

    # 5-fold rotation about the axis through raw[0] = (0, 1, phi)/|.|.
    axis5 = raw[0]
    R5 = _rotation_about_axis(axis5, 2.0 * np.pi / 5.0)

    # 3-fold rotation about a face center (choose center of face
    # containing raw[0], raw[2], raw[8]: average and normalise).
    face = (raw[0] + raw[2] + raw[8]) / 3.0
    axis3 = face / np.linalg.norm(face)
    R3 = _rotation_about_axis(axis3, 2.0 * np.pi / 3.0)

    # The group generated by the 5-fold, the 3-fold and inversion IS I_h —
    # 120 elements. The hand-rolled BFS that used to live here (with its own
    # dedup-by-allclose and its own iteration cap) was a second closure
    # implementation; the core's is the one that carries the cap as a NAMED
    # refusal.
    return list(_close_rigid_motions([R5, R3, _inversion_op()]))


def _rotation_about_axis(axis: np.ndarray, theta: float) -> RigidMotion:
    r"""Rotation about ``axis`` by ``theta``.

    Delegates to the core, which builds the rotation **plane** orthogonal to
    ``axis`` and rotates in it. "About an axis" is a three-dimensional
    convenience — the :math:`(d-2)`-dimensional fixed set happens to be a line
    — not a separate construction.
    """
    return RigidMotion.rotation_about_axis(axis=axis, angle=theta)


# ---------------------------------------------------------------------------
# Public attribute install
# ---------------------------------------------------------------------------

# Pre-instantiated singletons attached to the public class. These are
# the canonical way users obtain a SubgroupOfO3 for a named entry.
SubgroupOfO3.Trivial = SubgroupOfO3._from_named(_NamedSubgroup.Trivial)
SubgroupOfO3.Dinfh = SubgroupOfO3._from_named(_NamedSubgroup.Dinfh)
SubgroupOfO3.OctahedralOh = SubgroupOfO3._from_named(_NamedSubgroup.OctahedralOh)
SubgroupOfO3.IcosahedralIh = SubgroupOfO3._from_named(_NamedSubgroup.IcosahedralIh)
SubgroupOfO3.SO3 = SubgroupOfO3._from_named(_NamedSubgroup.SO3)
SubgroupOfO3.O3 = SubgroupOfO3._from_named(_NamedSubgroup.O3)


__all__ = [
    "AXIS_INDEX",
    "AXIS_LETTER",
    "IdentityComponent",
    "Realization",
    "SubgroupOfO3",
]
