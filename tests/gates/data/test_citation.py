r"""``Citation`` — the laws of a citation value (P1 step 2c, #405).

A :class:`~orpheus.data.citation.Citation` names a published work by its
``docs/refs.bib`` key and, optionally, a place in it. It is a VALUE: two
citations of the same place in the same work are one citation, so it must
behave as a dict key and a set member, it must be immutable, and it must be
parsed at construction (a key outside the BibTeX key grammar, or a blank
locator, is refused there and never reaches a registry).

Whether a key EXISTS in ``docs/refs.bib`` is not a law of the type (``docs/``
does not ship): it is the registry gate
``tests/gates/derivations/test_registry_citations_resolve.py``, which rests on
the grammar rows below.

Each refusal row pins its message fragment, and the two fragments are
disjoint (a blank-locator input must not read as a bad key, and the reverse),
so a refusal raised for the wrong reason is red. The refusal rows rest on the
positive row: a type that refused every input would pass them all.
"""

from __future__ import annotations

import dataclasses

import pytest

from orpheus.data.citation import Citation

pytestmark = pytest.mark.foundation

_SELF = "tests/gates/data/test_citation.py"

#: The refusal messages' distinctive fragments. Disjoint by design.
_BAD_KEY = "not a refs.bib key"
_BLANK_LOCATOR = "is blank"


def _require(condition: bool, message: str) -> None:
    """A ``-O``-surviving assertion (not a bare ``assert``)."""
    if not condition:
        pytest.fail(message)


# ── the positive leg ───────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "bibkey",
    ["SoodForsterParsons2003", "SoodLA13511_1999", "Atalay1997", "A", "a_1"],
)
def test_legal_bibkeys_construct(bibkey: str) -> None:
    """A key in the grammar (a letter, then letters, digits or underscores)
    constructs, with and without a locator, and keeps both fields verbatim."""
    whole = Citation(bibkey)
    placed = Citation(bibkey, "problem 4")
    _require(whole.bibkey == bibkey and whole.locator is None, f"{whole!r}")
    _require(placed.bibkey == bibkey and placed.locator == "problem 4", f"{placed!r}")


def test_citation_is_reexported_from_orpheus_data() -> None:
    """``orpheus.data.Citation`` is the one class, not a copy."""
    import orpheus.data

    _require(orpheus.data.Citation is Citation, "orpheus.data.Citation is another object")


# ── value semantics ────────────────────────────────────────────────────────


@pytest.mark.rests_on(f"{_SELF}::test_legal_bibkeys_construct")
@pytest.mark.parametrize(
    ("bibkey", "locator"),
    [("SoodForsterParsons2003", "problem 4"), ("Atalay1997", None)],
)
def test_equal_citations_are_one_value(bibkey: str, locator: str | None) -> None:
    """Two separately built citations of one place are equal, hash equal,
    collapse to one set member, and find each other's dict entry."""
    a, b = Citation(bibkey, locator), Citation(bibkey, locator)
    _require(a is not b, "the fixture must build two objects")
    _require(a == b, f"{a!r} != {b!r}")
    _require(hash(a) == hash(b), f"hash({a!r}) != hash({b!r})")
    _require(len({a, b}) == 1, f"the set {{a, b}} holds {len({a, b})} members")
    _require({a: "value"}[b] == "value", "a dict keyed by a does not find b")


@pytest.mark.rests_on(f"{_SELF}::test_legal_bibkeys_construct")
@pytest.mark.parametrize(
    ("fields_a", "fields_b"),
    [
        (("SoodForsterParsons2003", "problem 4"), ("SoodLA13511_1999", "problem 4")),
        (("SoodForsterParsons2003", "problem 4"), ("SoodForsterParsons2003", "problem 3")),
        (("SoodForsterParsons2003", None), ("SoodForsterParsons2003", "problem 4")),
    ],
    ids=["key-differs", "locator-differs", "whole-work-vs-place"],
)
def test_citations_differing_in_one_field_are_two_values(
    fields_a: tuple[str, str | None], fields_b: tuple[str, str | None]
) -> None:
    """A citation differing in one field is another value. Separation is
    asserted through the containers, never as ``hash(a) != hash(b)`` (equal
    hashes of unequal values are legal). Built in the body, never in the
    parametrize list, so a mutation that makes construction raise reddens
    rows instead of killing collection."""
    a, b = Citation(*fields_a), Citation(*fields_b)
    _require(a != b, f"{a!r} == {b!r}")
    _require(len({a, b}) == 2, f"the set {{a, b}} holds {len({a, b})} members")
    _require(b not in {a: "value"}, f"a dict keyed by {a!r} finds {b!r}")


@pytest.mark.rests_on(f"{_SELF}::test_legal_bibkeys_construct")
def test_locator_defaults_to_none() -> None:
    """``locator`` defaults to ``None``: the citation of the work as a whole."""
    whole = Citation("Atalay1997")
    _require(whole.locator is None, f"default locator is {whole.locator!r}")
    _require(whole == Citation("Atalay1997", None), "the default is not None")


@pytest.mark.rests_on(f"{_SELF}::test_legal_bibkeys_construct")
@pytest.mark.parametrize("field", ["bibkey", "locator"])
def test_citation_is_frozen(field: str) -> None:
    """Neither field can be assigned after construction."""
    citation = Citation("SoodForsterParsons2003", "problem 4")
    with pytest.raises(dataclasses.FrozenInstanceError):
        setattr(citation, field, "Atalay1997")
    _require(
        citation == Citation("SoodForsterParsons2003", "problem 4"),
        f"the refused assignment changed the value: {citation!r}",
    )


# ── the refusals, one per law ──────────────────────────────────────────────


@pytest.mark.rests_on(f"{_SELF}::test_legal_bibkeys_construct")
def test_refuses_empty_bibkey() -> None:
    """An empty key cites nothing."""
    with pytest.raises(ValueError, match=_BAD_KEY) as refused:
        Citation("", "problem 4")
    _require(_BLANK_LOCATOR not in str(refused.value), f"reads as a locator refusal: {refused.value}")


@pytest.mark.rests_on(f"{_SELF}::test_legal_bibkeys_construct")
@pytest.mark.parametrize(
    "bibkey",
    ["Sood Forster", "2003Sood", "_Sood2003", "Sood-2003", " Sood2003", "Sood2003\n"],
    ids=["space", "leading-digit", "leading-underscore", "hyphen", "leading-space", "trailing-newline"],
)
def test_refuses_illegal_bibkey(bibkey: str) -> None:
    """A key outside the grammar is refused. ``Sood-2003`` is the retired
    ``Provenance.paper_id`` spelling, which is not a key; the trailing newline
    reddens a ``re.match(... $)`` in place of ``fullmatch``."""
    with pytest.raises(ValueError, match=_BAD_KEY) as refused:
        Citation(bibkey)
    _require(_BLANK_LOCATOR not in str(refused.value), f"reads as a locator refusal: {refused.value}")


@pytest.mark.rests_on(f"{_SELF}::test_legal_bibkeys_construct")
@pytest.mark.parametrize("locator", ["", "   ", "\t\n"], ids=["empty", "spaces", "tab-newline"])
def test_refuses_blank_locator(locator: str) -> None:
    """A locator that is present must name a place; ``None`` is the spelling
    of "the whole work". The replace leg shows the law re-fires on
    ``dataclasses.replace`` (Pattern 4 corollary)."""
    with pytest.raises(ValueError, match=_BLANK_LOCATOR) as refused:
        Citation("SoodForsterParsons2003", locator)
    _require(_BAD_KEY not in str(refused.value), f"reads as a key refusal: {refused.value}")
    with pytest.raises(ValueError, match=_BLANK_LOCATOR):
        dataclasses.replace(Citation("SoodForsterParsons2003", "problem 4"), locator=locator)
