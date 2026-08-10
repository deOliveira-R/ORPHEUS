r"""The Peierls reference-naming rule obeys its own defining laws (#345).

:func:`~orpheus.derivations.continuous.peierls_nystrom.naming.reference_name`
is a *keying* function: it turns a reference's identity into the string every
consumer uses to find it — the registry key, the lazy-builder key, the
capability-matrix row label, and the ``ContinuousReferenceSolution.name``
stamp. Before #345 that rule was written out nine times across five modules,
and two of the spellings had already diverged.

A keying function has intrinsic laws, and a project rule that says every
math-bearing primitive ships a test of them rather than only of its usage.
The laws here:

**Injectivity.** Two distinct identities must not key to one name — otherwise
the registry silently loses a reference to a dict collision.

**Refusal of the unit error it exists to prevent.** ``r0_over_R`` is a
*dimensionless ratio*. The defect #345 removed was an absolute length flowing
into that slot, right only at :math:`R = 1`; the constructor now refuses
anything outside :math:`(0, 1)`, so the most common way to make that mistake
raises instead of producing a plausible-but-wrong name.

**Totality over the declared shapes, and refusal outside them.**

**Agreement with the tree.** The whole point of one rule is that the names it
produces are the names that exist, so the shipped grid is checked against the
live registry — cheaply, by name (materialising a Peierls reference is an
O(minutes) mpmath solve, Issue #212).
"""

from __future__ import annotations

import pytest

from orpheus.derivations.continuous.peierls_nystrom.cases import SHIPPED_CLASS_A
from orpheus.derivations.continuous.peierls_nystrom.naming import (
    SHAPE_TAGS,
    ShippedReference,
    reference_name,
)
from orpheus.derivations.reference_values import continuous_all_names

pytestmark = pytest.mark.foundation


class TestTheNameIsAFunctionOfIdentityAlone:
    """The keying laws, independent of what the tree happens to ship."""

    def test_it_is_injective_over_a_grid_that_varies_every_field(self) -> None:
        """Distinct identities key to distinct names.

        A collision does not raise anywhere — it silently overwrites a
        registry entry, and the lost reference is discoverable only by
        noticing the count is short.
        """
        identities = [
            ShippedReference(shape, n_groups=ng, n_regions=nr, r0_over_R=r0)
            for shape in SHAPE_TAGS
            for ng in (1, 2, 4)
            for nr in (1, 2, 4)
            for r0 in (None, 0.1, 0.2, 0.3, 0.55)
            # slab has no cavity; the ratio is not part of its identity
            if not (shape == "slab" and r0 is not None)
        ]
        names = [case.name for case in identities]
        collisions = {n for n in names if names.count(n) > 1}
        assert not collisions, (
            f"reference_name is not injective: {sorted(collisions)} is keyed "
            f"by more than one identity, so the registry would drop one."
        )
        assert len(identities) >= 30, f"grid too thin to be evidence: {len(identities)}"

    @pytest.mark.parametrize("bad", [1.0, 2.0, 0.0, -0.1, 100.0])
    def test_it_refuses_a_length_where_the_dimensionless_ratio_belongs(
        self, bad: float
    ) -> None:
        """The unit error #345 removed, made loud.

        ``1.0`` and ``2.0`` are the readings an absolute :math:`r_0`
        produces at a sub-unit outer radius; ``0.0`` is a solid shape, which
        is a Class-B case and a different builder.
        """
        with pytest.raises(ValueError, match="RATIO"):
            reference_name("cylinder-1d", n_groups=1, n_regions=1, r0_over_R=bad)

    def test_a_ratio_inside_the_shell_is_accepted(self) -> None:
        """The positive leg — refusal must not be the only behaviour tested.

        Without it the guard could reject everything and every negative
        test above would still pass (`vv-principles` #11).
        """
        assert reference_name(
            "sphere-1d", n_groups=2, n_regions=1, r0_over_R=0.25
        ) == "peierls_sph1D_hollow_2eg_1rg_r0_25"

    def test_an_unknown_shape_is_refused_by_name(self) -> None:
        with pytest.raises(ValueError, match="unknown shape"):
            reference_name("torus-3d", n_groups=1, n_regions=1)

    def test_every_declared_shape_keys(self) -> None:
        """Totality: the tag table and the function agree on the vocabulary."""
        for shape in SHAPE_TAGS:
            assert reference_name(shape, n_groups=1, n_regions=1).startswith(
                f"peierls_{SHAPE_TAGS[shape]}_"
            )

    def test_the_tag_rounds_the_RATIO_not_a_radius(self) -> None:
        """Pins the formula itself, which is where the divergence lived.

        `[M]` 2026-08-09: the capability matrix used ``round(r0 * 100)``
        while the registry used ``round(r0/R_out * 100)``. Both read ``10``
        at ``r0 = 0.1`` because every shipped ``R_out`` is ``1.0``. This
        asserts the surviving rule on a ratio no radius would produce.
        """
        assert reference_name(
            "cylinder-1d", n_groups=1, n_regions=1, r0_over_R=0.125
        ).endswith("_r0_12"), "the tag is round(100 * r0_over_R), zero-padded"


class TestTheRatioLawIsEnforcedAtTheBOUNDARY:
    """The same law, twice, on purpose — and the layer is the point.

    :func:`reference_name` refuses a non-ratio because a wrong tag is a
    wrong registry key. But it is called by the *stamp*, at the very END
    of a build, so on its own it means a caller pays an O(minutes) mpmath
    solve on a geometry whose cavity sits outside its own shell and then
    receives an error about **naming**.

    So ``build_two_surface_case`` refuses first. This is not a Pattern-2
    duplicated invariant: the two guards answer different questions — the
    boundary one protects the *caller's* time and states the API contract,
    the naming one protects the *registry key* for every other producer of
    a name. Both are cheap; neither subsumes the other.
    """

    @pytest.mark.parametrize("bad", [1.0, 2.0, 0.0, -0.1])
    @pytest.mark.parametrize("shape", ["cylinder-1d", "sphere-1d"])
    def test_a_non_ratio_is_refused_before_any_solve(
        self, shape: str, bad: float
    ) -> None:
        """Fast by construction — if this test is ever slow, the guard moved."""
        from orpheus.derivations.continuous.peierls_nystrom.cases import (
            build_two_surface_case,
        )

        with pytest.raises(ValueError, match="RATIO"):
            build_two_surface_case(shape, "1g", 1, r0_over_R=bad)

    @pytest.mark.parametrize("shape", ["cylinder-1d", "sphere-1d"])
    def test_a_missing_cavity_names_the_solid_alternative(self, shape: str) -> None:
        from orpheus.derivations.continuous.peierls_nystrom.cases import (
            build_two_surface_case,
        )

        with pytest.raises(ValueError, match="build_one_surface_compact_case"):
            build_two_surface_case(shape, "1g", 1)


class TestTheShippedGridIsTheTree:
    """One rule is only worth having if it names things that exist."""

    def test_every_grid_entry_is_a_registered_reference(self) -> None:
        """Cheap by name — ``continuous_all_names`` forces no builder.

        This is the identity half of the cross-check
        ``peierls_nystrom.rst`` promised. The *field* half stays ungated;
        see the module docstring of
        ``tests/derivations/test_capability_matrices.py`` for why.
        """
        registered = set(continuous_all_names())
        missing = sorted(c.name for c in SHIPPED_CLASS_A if c.name not in registered)
        assert not missing, (
            f"{len(missing)} shipped grid entries are not registered: {missing}\n"
            f"Either the registry walker stopped finding "
            f"continuous_case_builders(), or a grid entry was added without a "
            f"builder that can produce it."
        )

    def test_the_grid_has_no_duplicate_identities(self) -> None:
        assert len(set(SHIPPED_CLASS_A)) == len(SHIPPED_CLASS_A), (
            "SHIPPED_CLASS_A contains a repeated identity — the second one "
            "silently overwrites the first in continuous_case_builders()."
        )

    def test_the_grid_covers_what_it_claims(self) -> None:
        """Guards against a grid that quietly shrinks to nothing.

        A comprehension over an empty source makes every gate above pass
        vacuously (`vv-principles` #8, the signature-tautological class).

        ⚠ This is the *local* anti-vacuity guard for this module, not the
        primary pin on the inventory. That is
        ``test_continuous_registry_lazy.py::test_builder_keyset_is_the_shipped_class_a_inventory``,
        which asserts the full 13-name set against a **hand-written
        literal** authored independently of both the grid and
        :func:`reference_name` — strictly stronger, and the one to keep if
        anyone ever reads these two as duplicates.
        """
        assert len(SHIPPED_CLASS_A) == 13, (
            f"expected the 13 shipped Class-A references (1 slab + 2 shapes "
            f"× 3 cavity ratios × 2 group counts), found "
            f"{len(SHIPPED_CLASS_A)}. If a reference was genuinely added or "
            f"retired, update this count AND the capability matrix."
        )
        assert {c.shape for c in SHIPPED_CLASS_A} == set(SHAPE_TAGS)
