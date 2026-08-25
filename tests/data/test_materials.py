r"""``Materials`` — the stage-1 declaration (un-weld arc, charter R20/R21).

The declaration's own laws: admission (refuse-empty ONLY — a
cross-material refusal at admission would let a declared-but-unassigned
spectator flip admissibility, violating the leak principle), the
read-only mapping surface, identity semantics (``eq=False``), and
``restrict`` — the reachable-subset constructor that IS the
assigned-but-undeclared guard (its MaterialMesh discharge is pinned by
``tests/transport/test_material_mesh.py`` / the PR-TYPED-0 suite; here
the guard is pinned at its one home).

Positive AND negative per vv-principles #11 throughout.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.data.materials import Materials
from orpheus.derivations.common.xs_library import get_mixture

pytestmark = pytest.mark.foundation


def _require(condition: bool, message: str) -> None:
    """A ``-O``-firing assertion (NOT a bare ``assert``)."""
    if not condition:
        pytest.fail(message)


def _two() -> Materials:
    return Materials({0: get_mixture("A", "2g"), 3: get_mixture("B", "2g")})


class TestAdmission:
    def test_empty_declaration_refuses(self) -> None:
        with pytest.raises(ValueError, match="non-empty materials"):
            Materials({})

    def test_admission_is_assignment_independent(self) -> None:
        """A spectator with different content is LEGAL at declaration —
        coherence questions belong to the reads/mint over the set a
        caller restricts to (leak principle), never to admission."""
        a = get_mixture("A", "2g")
        b = get_mixture("A", "1g")  # different ng — still declarable
        _require(len(Materials({0: a, 1: b})) == 2, "superset declaration must admit")

    def test_mapping_is_read_only_after_admission(self) -> None:
        src = {0: get_mixture("A", "2g")}
        mats = Materials(src)
        with pytest.raises(TypeError):
            mats.mixtures[1] = get_mixture("B", "2g")  # type: ignore[index]
        # and the defensive copy means mutating the SOURCE dict is inert
        src[1] = get_mixture("B", "2g")
        _require(len(mats) == 1, "the declaration must not alias the source dict")


class TestMappingSurface:
    def test_reads(self) -> None:
        mats = _two()
        _require(mats[0].ng == 2, "__getitem__")
        _require(3 in mats and 1 not in mats, "__contains__")
        _require(sorted(mats) == [0, 3], "__iter__ over ids")
        _require(len(mats) == 2, "__len__")
        _require(mats.get(9) is None, ".get default")
        _require(set(mats.keys()) == {0, 3}, ".keys")
        _require(len(list(mats.values())) == 2, ".values")
        _require(dict(mats.items())[3] is mats[3], ".items")
        _require(mats.ids == frozenset({0, 3}), ".ids")

    def test_identity_semantics(self) -> None:
        """eq=False: two content-equal declarations are DISTINCT
        declarations (content identity joins the typed-axis identity
        family later)."""
        a, b = _two(), _two()
        _require(a == a, "reflexive identity")
        _require(a != b, "content-equal twins are distinct declarations")

    def test_of_is_parse_at_boundary(self) -> None:
        mats = _two()
        _require(Materials.of(mats) is mats, "an admitted declaration passes through")
        coerced = Materials.of({0: get_mixture("A", "2g")})
        _require(isinstance(coerced, Materials), "a bare dict is admitted once")


class TestRestrict:
    def test_subset_preserves_given_order_and_entries(self) -> None:
        mats = _two()
        sub = mats.restrict([3, 0])
        _require(list(sub) == [3, 0], "restrict preserves the given id order")
        _require(sub[3] is mats[3] and sub[0] is mats[0], "entries are the same objects")

    def test_spectator_is_silently_not_selected(self) -> None:
        """Declared-but-unassigned is legal and inert — no warning
        machinery exists, by ruling."""
        sub = _two().restrict([0])
        _require(sub.ids == frozenset({0}), "only the reachable id survives")

    def test_assigned_undeclared_refuses_naming_both_sets(self) -> None:
        with pytest.raises(ValueError, match=r"references material ids \[7\]"):
            _two().restrict([0, 7])

    def test_numpy_integer_ids_admit(self) -> None:
        """The pullback feeds np.unique output — numpy ints must key."""
        sub = _two().restrict(np.array([0, 3]))
        _require(sub.ids == frozenset({0, 3}), "np integer ids restrict cleanly")
