r""":class:`Descent` — the two realizations of an orbit space's functions, and the WITNESS that they are one.

The functions on :math:`M/H` have two honest spellings: **upstairs**, the
:math:`H`-invariant subspace of a basis on the base (for :math:`SO(2)_x` and
the real harmonics, the :math:`m = 0` column); **downstairs**, the quotient's
own classical family (:math:`\{P_\ell(\mu)\}`). They are isomorphic, and
without a witness they are a Cardinal-Rule-2 twin — two spellings of one
concept, free to drift. This module is the witness, and it is the
discriminator's gate: *which* spelling a frame binds is DERIVED from the
entry, never chosen at a call site.

Claim layer / pillar
--------------------

**Structural (L0)**, one named mutation each — not a value claim about any
solver. The one numerical row (:func:`test_the_two_realizations_agree_at_the_BIT`)
is a **regression contract** in the ``vv`` §bit-identity sense: the two sides
share ``scipy.special.lpmv``, which is legal below the trusted-library line
(``algebra-of-record``) but means this row on its own would test *convention
agreement*, not correctness. Its independent leg is
``tests/numerics/test_legendre_basis.py::test_legendre_values_against_an_independent_three_term_recurrence``
— read the two together (``vv-principles`` #22).

⚠ **The refusal is AXIS-keyed, and that is a measured design constraint.**
``[M]`` 2026-09-02, live descending slots per degree: about **x** it is 1 at
every :math:`\ell`; about **y** and **z** it is 1 at :math:`\ell \le 1` and
**0** from :math:`\ell \ge 2` (the invariant subspace is still
one-dimensional per degree by Schur — it simply stops being a set of SLOTS of
the harmonic table, since the harmonics' polar axis is x). So the three axes
are indistinguishable by slot structure at :math:`\ell \le 1`, and a refusal
keyed on *measured alignment* would be silently inert at exactly the order
where ``fission`` and ``n2n`` mint on every solve (:math:`L = 0`).
:func:`test_upstairs_face_is_refused_for_a_non_polar_axis` carries the
:math:`L = 1` row for that reason, and
:func:`test_the_measured_alignment_an_axis_keyed_refusal_must_ignore` records
the table.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.basis.descent import Descent
from orpheus.numerics.basis.indicator_basis import IndicatorBasis
from orpheus.numerics.basis.legendre_basis import LegendreBasis
from orpheus.numerics.basis.spherical_harmonic_basis import (
    MirrorEvenSphericalHarmonicBasis,
    SphericalHarmonicBasis,
)
from orpheus.numerics.manifold import ENERGY, SPHERE
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.symmetry import SubgroupOfO3

pytestmark = [pytest.mark.foundation]


_SPHERE_RULES: tuple[str, ...] = (
    "level_symmetric(4)",
    "level_symmetric(8)",
    "lebedev(11)",
    "lebedev(17)",
    "product(4,6)",
    "product(8,8)",
    "folded_product(4,8)",
)


def _rule(label: str) -> Quadrature:
    """A shipped rule from its LABEL — parametrize carries labels, never objects."""
    match label:
        case "level_symmetric(4)":
            return Quadrature.level_symmetric(4)
        case "level_symmetric(8)":
            return Quadrature.level_symmetric(8)
        case "lebedev(11)":
            return Quadrature.lebedev(11)
        case "lebedev(17)":
            return Quadrature.lebedev(17)
        case "product(4,6)":
            return Quadrature.product(n_mu=4, n_phi=6)
        case "product(8,8)":
            return Quadrature.product(n_mu=8, n_phi=8)
        case "folded_product(4,8)":
            return Quadrature.folded_product(n_mu=4, n_phi=8)
    if label.startswith("gauss_legendre("):
        return Quadrature.gauss_legendre(int(label[15:-1]))
    raise ValueError(f"unknown rule label {label!r}")


def _so2(axis: str):
    return SPHERE.quotient(SubgroupOfO3.SO2(axis))


# ══════════════════════════════════════════════════════════════════════
# The isomorphism, at the BIT
# ══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("label", _SPHERE_RULES)
def test_the_two_realizations_agree_at_the_BIT(label: str) -> None:
    r"""``downstairs(π(Ω)) == upstairs_columns(Ω)`` exactly — ``max|Δ| = 0.0``, 7 of 7 sphere rules at :math:`L = 4`.

    Stated at the BIT rather than "to machine precision" on purpose: it is a
    stronger and measured claim, and it is what lets #429's repair be a
    bit-identity regression contract on the 1-D flux instead of a tolerance.

    ``[M]`` 2026-09-02 at this tree, all seven rules, ``max_abs_difference``
    returns exactly ``0.0``.
    """
    nodes = np.asarray(_rule(label).measure.nodes, dtype=float)
    descent = Descent.for_entry(_so2("x"), 4)

    assert descent.max_abs_difference(nodes) == 0.0, (
        f"{label}: {descent.max_abs_difference(nodes):.3e}"
    )
    assert descent.is_isomorphism(nodes)
    # the shape is the honest denominator: ONE column per degree, not one
    # per table slot (the |m| > l padding descends vacuously — B5b).
    assert descent.upstairs_columns(nodes).shape == (nodes.shape[0], 5)


@pytest.mark.parametrize("label", ("gauss_legendre(8)", "gauss_legendre(16)"))
@pytest.mark.parametrize("L", (0, 1, 2))
def test_the_isomorphism_holds_on_the_one_d_rules_lifted_to_the_base(
    label: str, L: int
) -> None:
    r"""The same identity on the 1-D rules, at :math:`L \le 2`.

    A 1-D rule's nodes are cosines, so to ask the UPSTAIRS side anything at
    all they must be lifted to the base. The only honest lift is a
    representative on :math:`S^2`; this row builds one that IS on the sphere
    (:math:`(\mu, \sqrt{1-\mu^2}, 0)`), which is what makes the comparison
    legitimate — the forged :math:`(\mu, 0, 0)` of ERR-080 is NOT on the
    sphere and 0.6 refuses it.

    ⭐ That the identity survives a lift whose azimuth is **not** zero is the
    content: :math:`P_\ell(\mu)` does not depend on the representative, which
    is precisely the property ERR-080's construction violated.
    """
    mu = np.asarray(_rule(label).measure.nodes, dtype=float).reshape(-1)
    lifted = np.column_stack([mu, np.sqrt(np.maximum(0.0, 1.0 - mu**2)), np.zeros_like(mu)])

    descent = Descent.for_entry(_so2("x"), L)
    assert descent.max_abs_difference(lifted) == 0.0, f"{label} L={L}"

    # …and the pullback really did go through pi: the downstairs values are
    # the same for a DIFFERENT representative of the same orbit.
    other = np.column_stack([mu, np.zeros_like(mu), np.sqrt(np.maximum(0.0, 1.0 - mu**2))])
    down = descent.downstairs
    assert down is not None
    assert np.array_equal(
        down.evaluate(_so2("x").quotient_map(lifted)),
        down.evaluate(_so2("x").quotient_map(other)),
    )


# ══════════════════════════════════════════════════════════════════════
# B7 — the axis-keyed refusal, and its L <= 1 row
# ══════════════════════════════════════════════════════════════════════


def test_upstairs_face_is_refused_for_a_non_polar_axis() -> None:
    r"""B7 — ``SO2("x")`` is admitted upstairs; ``SO2("y")`` / ``("z")`` are refused, naming the work.

    ⛔ **The refusal must be keyed on the AXIS, never on measured slot
    alignment.** ``[M]`` 2026-09-02: about y and z the invariant subspace IS
    slot-aligned at :math:`\ell \le 1` (z, :math:`\ell=1` → column 0; y,
    :math:`\ell=1` → column 2) and only spreads from :math:`\ell \ge 2`. An
    alignment-keyed refusal would therefore be **silently inert at**
    :math:`L \le 1` — the order at which ``fission.py`` and ``n2n.py`` mint a
    moment space on *every* solve, isotropic included. So the row below
    exercises :math:`L = 1` as well as :math:`L = 3`: at :math:`L = 1` an
    alignment-keyed implementation passes and this gate must still red it.
    """
    nodes = np.asarray(_rule("level_symmetric(4)").measure.nodes, dtype=float)

    # positive leg — the harmonics' own polar axis
    for L in (1, 3):
        assert Descent.for_entry(_so2("x"), L).upstairs_columns(nodes).shape == (
            nodes.shape[0],
            L + 1,
        )

    # negative legs — including the L = 1 row an alignment key would admit
    for axis in ("y", "z"):
        for L in (1, 3):
            descent = Descent.for_entry(_so2(axis), L)
            with pytest.raises(NotImplementedError, match="polar axis"):
                descent.upstairs_columns(nodes)
            with pytest.raises(NotImplementedError, match="no masked harmonic basis"):
                _ = descent.upstairs
            # …and the DOWNSTAIRS face is available at every axis, which is
            # why the refusal costs nothing (it names work, not a wall).
            assert descent.downstairs == LegendreBasis(L=L, axis=axis)


def test_the_measured_alignment_an_axis_keyed_refusal_must_ignore() -> None:
    r"""The measurement that FORCES B7's key — recorded here so the design cannot be relaxed by inspection.

    About y and z the :math:`SO(2)`-invariant subspace is one-dimensional per
    degree (Schur), and this row measures how many SLOTS of the harmonic
    table carry it. ``[M]`` 2026-09-02, live descending slots per degree:

    ======  ===================  ===================
    axis    :math:`\ell \le 1`    :math:`\ell \ge 2`
    ======  ===================  ===================
    x       1                    1
    y, z    1                    **0**
    ======  ===================  ===================

    So at :math:`L \le 1` the three axes are INDISTINGUISHABLE by slot
    structure, and only from :math:`\ell \ge 2` does the y/z subspace stop
    being a set of slots at all (it spreads over several, none of them
    fibre-constant on its own). An implementation that refused on measured
    alignment would therefore admit y and z at :math:`L \le 1` — the order
    at which every solve mints — which is why B7's key is the axis.
    """
    harmonic_axis = _so2("x")
    for L in (1, 4):
        live = SphericalHarmonicBasis(L=L).live_slot_mask
        polar = (harmonic_axis.descending_slots(SphericalHarmonicBasis(L=L)) & live).sum(
            axis=1
        )
        assert polar.tolist() == [1] * (L + 1), (
            f"about x, L={L}: per-degree live descending counts {polar.tolist()}"
        )
        for axis in ("y", "z"):
            other = (
                _so2(axis).descending_slots(SphericalHarmonicBasis(L=L)) & live
            ).sum(axis=1)
            expected = [1] * min(L + 1, 2) + [0] * max(0, L - 1)
            assert other.tolist() == expected, (
                f"axis {axis}, L={L}: per-degree live descending counts "
                f"{other.tolist()}"
            )
            if L == 1:
                assert other.tolist() == polar.tolist(), (
                    "at L=1 the three axes agree — the measurement that makes an "
                    "alignment-keyed refusal inert"
                )
            else:
                assert other.tolist() != polar.tolist()


# ══════════════════════════════════════════════════════════════════════
# The DISCRIMINATOR — which realization a frame binds
# ══════════════════════════════════════════════════════════════════════


def test_the_discriminator_picks_downstairs_iff_the_quotient_has_a_classical_basis() -> None:
    r"""The sentence, as a property: Legendre for :math:`SO(2)`, σ-even harmonics for a mirror, the parent for Trivial.

    This is the single-source claim the class exists to make: the basis a
    frame on an orbit space carries is DERIVED from the entry, so a call site
    cannot pick a different one and no second copy of the rule exists.
    """
    legendre_entry = Descent.for_entry(_so2("x"), 3)
    assert legendre_entry.downstairs == LegendreBasis(L=3, axis="x")
    assert legendre_entry.frame_basis == LegendreBasis(L=3, axis="x")

    mirror_entry = Descent.for_entry(SPHERE.quotient(SubgroupOfO3.Mirror("y")), 3)
    assert mirror_entry.downstairs is None, (
        "a mirror quotient has NO classical named family — that is what makes "
        "its upstairs face the only spellable realization"
    )
    assert mirror_entry.frame_basis == MirrorEvenSphericalHarmonicBasis(L=3, mirror_axis=1)

    trivial_entry = Descent.for_entry(SPHERE.quotient(SubgroupOfO3.Trivial), 3)
    assert trivial_entry.downstairs == SphericalHarmonicBasis(L=3)
    assert trivial_entry.frame_basis == SphericalHarmonicBasis(L=3)

    # the discriminator states itself, and the sentence names both branches
    sentence = legendre_entry.discriminator
    assert "downstairs" in sentence and "upstairs" in sentence

    # NEGATIVE: a mirror entry must NOT answer Legendre, and an SO(2) entry
    # must NOT answer the mirror basis — the two branches are exclusive.
    assert not isinstance(mirror_entry.frame_basis, LegendreBasis)
    assert not isinstance(
        legendre_entry.frame_basis, MirrorEvenSphericalHarmonicBasis
    )


def test_max_abs_difference_refuses_an_entry_with_no_downstairs_face() -> None:
    r"""A mirror entry has no downstairs realization, so there is nothing to compare — it says so."""
    nodes = np.asarray(_rule("level_symmetric(4)").measure.nodes, dtype=float)
    descent = Descent.for_entry(SPHERE.quotient(SubgroupOfO3.Mirror("y")), 2)
    with pytest.raises(NotImplementedError, match="no\n?\\s*downstairs realization|downstairs realization"):
        descent.max_abs_difference(nodes)


# ══════════════════════════════════════════════════════════════════════
# Construction refusals
# ══════════════════════════════════════════════════════════════════════


def test_for_entry_refuses_a_quotient_of_anything_but_the_sphere() -> None:
    r"""``for_entry`` supplies the real harmonics as the parent, so only sphere entries have one."""
    energy_entry = ENERGY.quotient(SubgroupOfO3.Trivial)
    with pytest.raises(NotImplementedError, match="only quotients of S\\^2"):
        Descent.for_entry(energy_entry, 2)

    # positive leg: a sphere entry constructs and carries the harmonics
    descent = Descent.for_entry(_so2("x"), 2)
    assert descent.parent == SphericalHarmonicBasis(L=2)


def test_a_parent_that_does_not_eat_the_base_is_refused_at_construction() -> None:
    r"""A descent pulls functions on the orbit space back to its BASE — a parent on another manifold is not one."""
    entry = _so2("x")
    foreign = IndicatorBasis(
        edges_per_axis=(np.array([0.0, 0.5, 1.0]),), partition_of=ENERGY
    )
    with pytest.raises(ValueError, match="a descent pulls functions"):
        Descent(entry=entry, parent=foreign)

    # positive leg
    Descent(entry=entry, parent=SphericalHarmonicBasis(L=2))


def test_a_parent_with_no_truncation_order_is_refused_when_the_order_is_needed() -> None:
    r"""``downstairs`` must know :math:`L`; a parent carrying no ``L`` is refused by TYPE, not by ``AttributeError``."""
    entry = _so2("x")
    edges = tuple(np.linspace(-1.0, 1.0, 4) for _ in range(3))
    parent = IndicatorBasis(edges_per_axis=edges, partition_of=SPHERE)
    descent = Descent(entry=entry, parent=parent)
    with pytest.raises(TypeError, match="carries no truncation order"):
        _ = descent.downstairs
