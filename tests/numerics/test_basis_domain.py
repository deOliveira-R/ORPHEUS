r"""Every basis states the manifold its functions EAT — #429 tracker 2.1.

A basis function is a **map**, and a map is not defined until its source is:
:math:`Y_\ell^m : S^2 \to \mathbb R` takes a POINT of :math:`S^2`;
:math:`\mathbf 1_R` takes a point of whatever was partitioned. Until 2.1 the
:class:`~orpheus.numerics.basis.base.Basis` ABC had no way to ask, so the
answer was smuggled through a coefficient-space NAME STRING — and one of the
two producers hard-coded it.

The defect that forced this
===========================

:meth:`~orpheus.numerics.basis.indicator_basis.IndicatorBasis.space` built
``f"L2[coarse_cells_R{ndim}]"``, asserting a SPATIAL manifold whatever the
edges actually partitioned. `[M]` 2026-09-01, pre-fix: a 2-group **energy**
indicator space and a 2-cell **spatial** indicator space were ``==``-equal
**and hash-equal**, both ``FunctionSpace('L2[coarse_cells_R1]', shape=(2,))``.
:class:`~orpheus.numerics.space.FunctionSpace` identity is ``(name, shape)``,
so a false name is not cosmetic — it is an illegal state that IS
representable, and every composability guard downstream reads it as truth.

⭐ **The sharpest form of it, and the reason this file's keystone is
``test_d6``:** at four of the five production sites the basis and its measure
are built in the SAME function, three to five lines apart, and the *measure*
named the manifold correctly the whole time (``support="energy"``,
``"spatial_R1"``, ``f"index({label})"``). The answer was never unavailable —
only unasked. So the durable gate is not "the name is right" but **"the two
halves of one frame name ONE manifold"**, which goes red the day either side
re-invents it.

What is NOT claimed here
========================

``DiscreteMeasure.support`` is still a ``str`` (that retype is tracker 2.0c),
so ``measure.py``'s ``f"L2[{support}]"`` and the basis's new
``f"L2[coarse_cells({domain.name})]"`` remain two spellings of the L² name.
``test_d6`` pins their AGREEMENT, which is the property that matters and the
thing 2.0c must preserve.

See :doc:`/theory/foundations/manifolds` §(b) and
:doc:`/theory/foundations/spaces`.
"""

from __future__ import annotations

import numpy as np
import pytest
from numpy.typing import NDArray

from orpheus.data.energy_grid import EnergyGrid
from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.basis import (
    Basis,
    IndicatorBasis,
    OverlapBasis,
    SphericalHarmonicBasis,
    WeightedIndicatorBasis,
)
from orpheus.numerics.basis.spherical_harmonic_basis import (
    MirrorEvenSphericalHarmonicBasis,
)
from orpheus.numerics.manifold import (
    ENERGY,
    SPHERE,
    EnergyGroups,
    IndexSet,
    Manifold,
    RealSpace,
    ambient_dim,
)
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.space import FunctionSpace
from orpheus.numerics.symmetry import SubgroupOfO3
from orpheus.sn.operators.loss_kernel_gauge import LossKernelBasis

pytestmark = [pytest.mark.foundation]


_TWO_CELLS = (np.array([0.0, 0.5, 1.0]),)
_TWO_GROUPS = (np.arange(3, dtype=float) - 0.5,)
_THREE_CELLS = (np.array([0.0, 0.3, 0.6, 1.0]),)


def test_d1_an_energy_space_and_a_spatial_space_are_not_the_same_space() -> None:
    """⭐⭐ KEYSTONE — the defect 2.1 exists to make unspellable.

    Two partitions with the same cell COUNT of two entirely different
    manifolds must not produce equal coefficient spaces. `[M]` before 2.1 they
    did, ``==`` and hash-equal alike.

    Both legs matter (``vv-principles`` #11). The negative control — a 3-cell
    spatial partition against the 2-cell one — is what shows the inequality is
    not a degenerate always-false: if ``__eq__`` were broken outright, the
    first assertion would pass for the wrong reason and this file would
    certify nothing.
    """
    energy = IndicatorBasis(_TWO_GROUPS, EnergyGroups(2)).space
    spatial = IndicatorBasis(_TWO_CELLS, RealSpace(1)).space

    assert energy.shape == spatial.shape == (2,)     # the collision's precondition
    assert energy != spatial
    assert hash(energy) != hash(spatial)

    # negative control: inequality is decided by the MANIFOLD, and shape still
    # separates two partitions of the SAME manifold.
    three = IndicatorBasis(_THREE_CELLS, RealSpace(1)).space
    assert spatial != three
    assert IndicatorBasis(_TWO_CELLS, RealSpace(1)).space == spatial


def test_d2_every_shipped_basis_answers_what_its_functions_eat() -> None:
    """The completeness claim, enumerated by RUNTIME introspection.

    ⚠ Deliberately **not** an AST census and not a hand-written list. `[M]`
    2026-09-01 an AST pass over this tree reported 3 direct / 5 recursive
    ``Basis`` subclasses where the runtime answer is **4 / 6** — inheritance is
    a runtime relation, and a static pass cannot see re-exports, aliased bases
    or qualified base names. Enumerating at runtime also means a basis added
    tomorrow is in scope without editing this file.

    Filtered to ``orpheus.``-defined classes so a test-local subclass declared
    by whichever module pytest imported first cannot make the gate
    order-dependent.
    """

    def shipped(cls: type) -> list[type]:
        out: list[type] = []
        for sub in cls.__subclasses__():
            if sub.__module__.startswith("orpheus."):
                out.append(sub)
            out += shipped(sub)
        return out

    subclasses = shipped(Basis)
    assert len(subclasses) >= 6, (
        f"expected at least the 6 shipped bases, found {len(subclasses)}: "
        f"{[c.__name__ for c in subclasses]} — if a basis was retired, retire "
        f"its row here too; if the count COLLAPSED, an import is missing and "
        f"this gate is measuring nothing."
    )
    for sub in subclasses:
        assert "domain" not in getattr(sub, "__abstractmethods__", frozenset()), (
            f"{sub.__name__} does not say what its functions eat"
        )


def test_d3_a_basis_that_cannot_say_what_it_eats_cannot_be_built() -> None:
    """The refusal is STRUCTURAL — at construction, not at first call.

    The stub implements every other abstract member, so ``domain`` is the only
    thing missing and the refusal cannot be credited to the wrong arm
    (``vv-principles`` #17's granularity trap). The positive leg is the same
    stub WITH ``domain``, which must construct.
    """

    class _Incomplete(Basis):
        # Full signatures, and bodies that RAISE rather than `...`: a stub
        # whose members silently return None is not a faithful stand-in for
        # the ABC, and the type checker says so.  Nothing here is ever
        # reached — the class cannot be instantiated, which is the claim.
        def evaluate(self, points: NDArray, /) -> NDArray:
            raise NotImplementedError

        def synthesize(self, c: NDArray, table: NDArray, /) -> NDArray:
            raise NotImplementedError

        def analyze(self, v: NDArray, table: NDArray, w: NDArray, /) -> NDArray:
            raise NotImplementedError

        def analyze_transpose(
            self, c: NDArray, table: NDArray, w: NDArray, /,
        ) -> NDArray:
            raise NotImplementedError

        def reconstruct(self, c: NDArray, table: NDArray, /) -> NDArray:
            raise NotImplementedError

        def reconstruct_transpose(
            self, v: NDArray, table: NDArray, /,
        ) -> NDArray:
            raise NotImplementedError

        def mass_matrix(self, measure: DiscreteMeasure, /) -> NDArray:
            raise NotImplementedError

        @property
        def space(self) -> FunctionSpace:
            raise NotImplementedError

    # The ignore DOCUMENTS the claim rather than hiding a defect: pyright
    # reports exactly the refusal this line asserts ("Basis.domain is not
    # implemented"), which is the established idiom for a negative leg here
    # (e.g. tests/transport/test_timed_full_field.py:126).
    with pytest.raises(TypeError, match="domain"):
        _Incomplete()                                # type: ignore[abstract]

    class _Complete(_Incomplete):
        @property
        def domain(self) -> Manifold:
            return SPHERE

    assert _Complete().domain == SPHERE                  # positive leg


def test_d4_a_partition_must_have_one_axis_per_coordinate() -> None:
    """The construction invariant, both legs.

    A ``d``-axis tensor partition partitions a ``d``-coordinate manifold. The
    negative leg is what stops a 2-axis spatial partition from claiming the
    energy axis; the positive legs are every shipped combination, so the guard
    cannot pass by refusing everything.

    ⚠ The invariant is the AMBIENT width, NOT
    :meth:`~orpheus.numerics.manifold.Manifold.contains` on the cell centres.
    The stronger check reads better and is wrong: `[M]` the single-region index
    partition ``[-0.5, n-0.5]`` that ``frame.py``'s axis marginal ships has
    centre :math:`(n-1)/2`, not an integer, and ``IndexSet`` admits only
    integers — so it would refuse a correct production caller
    (``vv-principles`` #16).
    """
    ok = [
        (_TWO_CELLS, RealSpace(1)),
        (_TWO_GROUPS, EnergyGroups(2)),
        (_TWO_GROUPS, IndexSet("g", 2)),
        ((np.array([0.0, 1.0]), np.array([0.0, 1.0, 2.0])), RealSpace(2)),
    ]
    for edges, manifold in ok:
        basis = IndicatorBasis(edges, manifold)
        assert ambient_dim(basis.domain) == basis.ndim

    with pytest.raises(ValueError, match="partition axis"):
        IndicatorBasis((np.array([0.0, 1.0]), np.array([0.0, 1.0])), ENERGY)
    with pytest.raises(ValueError, match="partition axis"):
        IndicatorBasis(_TWO_CELLS, RealSpace(2))
    with pytest.raises(ValueError, match="partition axis"):
        IndicatorBasis(_TWO_CELLS, SPHERE)               # 1 axis, 3 coordinates


@pytest.mark.parametrize("L", [0, 1, 3, 7])
def test_d5_a_harmonic_eats_a_direction_at_every_degree(L: int) -> None:
    """:math:`S^2` is constant in :math:`L` — truncation changes the SPAN, not
    the source.

    Swept over four degrees rather than asserted once, because "the domain does
    not depend on ``L``" is the actual claim and a single degree cannot carry
    it.
    """
    assert SphericalHarmonicBasis(L=L).domain == SPHERE


@pytest.mark.parametrize("axis,name", [(0, "S^2/sigma_x"), (1, "S^2/sigma_y"), (2, "S^2/sigma_z")])
def test_d5b_a_mirror_even_harmonic_eats_the_QUOTIENT(axis: int, name: str) -> None:
    """⭐ The class's whole content, stated as a type.

    Every σ-even harmonic takes the same value at :math:`\\Omega` and
    :math:`\\sigma_a\\Omega`, and a function constant on the orbits of
    :math:`H` is a function on :math:`M/H`. So this basis's domain is the
    quotient — and it must NOT compare equal to the sphere, which is the
    discrimination that had no operands before 2.1.

    `[M]` the derived name matches the support a ``folded_product`` rule's
    angular frame already reports (landed 0.1c), so the two halves of a folded
    frame agree by DERIVATION rather than by two tags that happen to match.
    """
    basis = MirrorEvenSphericalHarmonicBasis(L=2, mirror_axis=axis)
    expected = SPHERE.quotient(SubgroupOfO3.Mirror("xyz"[axis]))

    assert basis.domain == expected
    assert basis.domain.name == name
    assert basis.domain != SPHERE                    # the discriminating leg
    assert SphericalHarmonicBasis(L=2).domain == SPHERE   # the parent is unmoved


def test_d6_a_frames_two_halves_name_ONE_manifold() -> None:
    """⭐⭐ THE FLAGSHIP. The basis and the measure of one frame agree.

    This is the property whose absence WAS the defect: at every site below the
    two objects are built within five lines of each other, and before 2.1 the
    measure named the manifold correctly while the basis hard-coded a spatial
    one. A gate on the basis's name alone would have been satisfied by any
    self-consistent lie; this one cannot be, because the two names come from
    independently-authored halves.

    ⚠ Two known limits, stated rather than hidden:

    * the fifth production pair (``frame.py``'s axis marginal, a private
      ``_collapse_pair`` whose frame is deliberately forgetful) is pinned
      instead by ``tests/numerics/test_axis_marginal.py``'s independent
      re-spelling of the same construction;
    * ``LossKernelBasis`` agrees on the point SET and not on its spelling —
      its measure tags the bare label while ``IndexSet`` wraps it as
      ``index(...)``. That divergence is asserted below rather than skipped,
      so tracker **2.0c** (which retypes ``support`` to a ``Manifold``) must
      come back here and cannot resolve it by accident.
    """
    mesh = Mesh1D(
        edges=np.array([0.0, 0.5, 1.0]), mat_ids=np.array([0, 1]),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"), bc_right=BC("vacuum"),
    )
    grid = EnergyGrid(edges=np.array([2.0e7, 1.0, 1.0e-5]))

    pairs = [
        ("Mesh1D", mesh.indicator_basis(), mesh.volume_measure),
        ("EnergyGrid", grid.as_basis(), grid.as_measure()),
    ]
    for who, basis, measure in pairs:
        assert basis.domain.name == measure.support, (
            f"{who}: the frame's basis says it lives on "
            f"{basis.domain.name!r} and its measure says {measure.support!r} — "
            f"one frame, two manifolds."
        )

    # ...and the spaces those agreeing halves mint stay distinguishable.
    assert mesh.indicator_basis().space != grid.as_basis().space

    # The loss-kernel pair: ONE point set, TWO spellings — pinned, not skipped.
    block = LossKernelBasis(table=np.eye(2), orbit=(0, 1), group=3)
    assert isinstance(block.domain, IndexSet)
    assert block.domain.label == "sn_trace_orbit(0, 1)_g3"
    assert block.domain.name == "index(sn_trace_orbit(0, 1)_g3)"


def test_d7_a_wrapping_basis_cannot_drift_from_the_basis_it_wraps() -> None:
    """Delegation, not duplication (Cardinal Rule 2).

    ``WeightedIndicatorBasis`` reweights the trial indicator, and ``OverlapBasis``
    decorates it with a fractional table; neither moves the support, so neither
    may carry its own answer. Asserted by IDENTITY (``is``), which a copied
    literal would fail even when the copy happens to be right today.
    """
    trial = IndicatorBasis(_TWO_GROUPS, EnergyGroups(2))

    weighted = WeightedIndicatorBasis(trial, np.ones(2))
    assert weighted.domain is trial.domain
    assert weighted.space == trial.space

    table = np.array([[1.0, 0.0], [0.0, 1.0]])
    overlap = OverlapBasis.from_indicator(trial, table)
    assert overlap.domain is trial.domain
    assert overlap.partition_of is trial.partition_of
