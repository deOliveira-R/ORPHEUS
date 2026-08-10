"""Gate: no docstring names a Python object that does not exist.

An unresolved Python-domain cross-reference (``:func:`` / ``:class:`` /
``:meth:`` / ``:mod:`` / ``:attr:``) renders as **plain text** and draws no
``-W`` warning. The natural next move — a nitpicky (``-n``) build — does not
help either, and the reason is structural rather than a config gap:

**Sphinx can only diagnose what it RENDERS.** A docstring in a module that is
never ``automodule``'d is never read, so its roles never become nodes and
nothing can warn at any severity. This tree carries 45 live ``automodule``
directives against 263 modules in ``orpheus/``, and ``tests/`` is never
rendered at all while holding ~500 fully-qualified xrefs.

So this gate resolves each target by **importing** it, via
``tools/check_docstring_xrefs.py``. Render coverage is irrelevant.

⛔ **CORRECTED 2026-08-10.** This docstring previously claimed the project's
knowledge graph "inherits the same blind spot for the same reason … and it
cannot be gated on." **That is false, and it was false in the safe-looking
direction** — it excused not reconciling the two instruments. Measured: the
five ``orpheus/derivations/continuous/`` reserved placeholders are README-only
directories; Nexus called all five dead and was right, while this checker
called all five alive and was wrong (``importlib`` materialises a PEP 420
namespace package from a bare directory). The instruments have *complementary*
blind spots — see the checker's module docstring — and the set difference
between them is a triage queue, not noise.

What is gated here, in two layers
---------------------------------

**End-to-end**, per tree, in a subprocess. It must be a subprocess: a full scan
imports arbitrary modules while probing, and doing that inside the pytest
process would pollute ``sys.modules`` for every test that runs after it.

**Per stage**, in-process, on real tree objects. A green end-to-end run cannot
distinguish "every reference resolves" from "the scanner found nothing"
(``vv-principles`` #17 — a broken instrument fails OPEN and reads as *write
more tests*). These gates therefore drive :func:`judge` directly and include a
positive control: a deliberately-dead target that MUST be reported. The
in-process import footprint is a handful of named modules the suite already
imports, not the arbitrary set a full scan touches — that is why the
subprocess rule above applies to the scan and not to these.

When this gate reds, the fix is per-site adjudication, not a blanket delete —
see ``.claude/rules/coding-standards.md`` "Retire as you go": a past-tense
mention is correct history and becomes a ``literal``; a present-tense claim or
an imperative instruction is a MUST-FIX. Do NOT silence a hit by adding to the
checker's ``ALLOWLIST``, which ships empty deliberately.
"""

from __future__ import annotations

import pathlib
import re
import subprocess
import sys

import pytest

from tools.check_docstring_xrefs import (
    Outcome,
    _is_empty_namespace_package,
    _self_attributes,
    extract_target,
    iter_text_blocks,
    judge,
    module_name_of,
    resolve,
)

REPO_ROOT = pathlib.Path(__file__).parent.parent
CHECKER = REPO_ROOT / "tools" / "check_docstring_xrefs.py"

pytestmark = pytest.mark.foundation


@pytest.mark.parametrize("tree", ["orpheus", "tests", "docs"])
def test_no_docstring_names_a_missing_object(tree: str) -> None:
    """Every decidable xref in ``tree`` resolves against the live namespace.

    Parametrised by tree so a production regression, a test-tree regression and
    a corpus regression are distinguishable in the report rather than sharing
    one red row.

    ``docs`` is a weaker case than the other two — a page under ``docs/`` IS
    rendered, so a nitpicky build could in principle flag its dead roles. It is
    gated here anyway for one gate instead of two and no ``nitpick_ignore``
    curation; see the checker's module docstring.
    """
    result = subprocess.run(
        [sys.executable, str(CHECKER), tree],
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )
    assert result.returncode == 0, (
        f"Dead Sphinx cross-references in {tree}/ — each one names a Python "
        f"object that no longer exists. No Sphinx build can see these (see this "
        f"module's docstring for why), so this gate is the only thing standing "
        f"between a module move and a docstring that lies.\n\n"
        f"Adjudicate each site (repoint / past-tense literal / rewrite), do not "
        f"blanket-delete, and do not add to the checker's ALLOWLIST.\n\n"
        f"{result.stdout}{result.stderr}"
    )


def test_the_checker_still_decides_most_of_what_it_finds() -> None:
    """A coverage FLOOR, so a silent collapse back to absolute-only reds.

    The checker judged 48.9 % of the roles it found until relative resolution
    landed (2026-08-10) and 71.7 % after. Nothing else in the suite would
    notice if that regressed: coverage falls, every gate above stays green, and
    the tree reads clean because half the surface stopped being looked at.

    The floor is deliberately well below the measured value — this catches a
    structural collapse, not ordinary drift from new prose.
    """
    result = subprocess.run(
        [sys.executable, str(CHECKER), "--quiet"],
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )
    def summary_line(label: str) -> int:
        match = re.search(rf"{label}\s*:\s*(\d+)", result.stdout)
        assert match is not None, (
            f"the checker's summary no longer reports {label!r} — this gate reads "
            f"its stdout, so a report-format change silently disables it.\n"
            f"{result.stdout}{result.stderr}"
        )
        return int(match.group(1))

    found = summary_line("xref roles found")
    decided = summary_line(r"decidable \(gated\)")
    assert found > 10_000, f"the scanner found only {found} roles — is it reading?"
    assert decided / found > 0.60, (
        f"the checker now decides only {decided}/{found} = {decided/found:.1%} of the "
        f"roles it finds, down from a measured 71.7 %. Relative-namespace "
        f"resolution has probably regressed (see `judge` / `candidate_paths`), "
        f"which silently un-gates thousands of references."
    )


class TestTheInstrumentHasTeeth:
    """The positive controls. A gate that cannot red is not a gate.

    Each leg names a target whose status is decided by reading the tree, so
    these do not drift with the corpus.
    """

    def test_a_dead_relative_target_is_REPORTED(self) -> None:
        """The control: a method that does not exist, in a real namespace.

        If this ever passes as ALIVE or DECLINED, the relative-resolution path
        is dead and the +2965 roles it gates are silently unchecked.
        """
        verdict = judge(
            "EigenvalueSolver.no_such_method_exists_here",
            ("orpheus.numerics.eigenvalue",),
        )
        assert verdict.outcome is Outcome.DEAD, verdict
        assert verdict.spelling == (
            "orpheus.numerics.eigenvalue.EigenvalueSolver.no_such_method_exists_here"
        )

    def test_a_live_relative_target_is_ALIVE(self) -> None:
        """The positive leg — refusal must not be the only behaviour tested.

        Without it the resolver could report everything dead and the control
        above would still pass (``vv-principles`` #11).
        """
        verdict = judge(
            "EigenvalueSolver.compute_keff", ("orpheus.numerics.eigenvalue",)
        )
        assert verdict.outcome is Outcome.ALIVE, verdict

    def test_a_target_imported_from_elsewhere_is_ALIVE(self) -> None:
        """Why relative resolution works at all: ``import`` IS the declaration.

        ``face_layout`` does not define ``Hashable``; it imports it. So the name
        is bound on the module and ``getattr`` finds it — which is exactly why
        this tool can read relative prose without emulating Sphinx.
        """
        verdict = judge("Hashable", ("orpheus.numerics.face_layout",))
        assert verdict.outcome is Outcome.ALIVE, verdict

    @pytest.mark.parametrize(
        "target",
        [
            "NotImplementedError",  # builtin, single-segment
            "TypeError",
            "dict.get",  # builtin, with an attribute walk
            "mpmath.quad",  # 3rd-party — 28 sites the curated root list missed
            "functools.lru_cache",  # stdlib
            "typing.overload",
        ],
    )
    def test_a_builtin_or_foreign_target_is_DECIDED(self, target: str) -> None:
        """⛔ This leg previously asserted these were DECLINED. That was the bug.

        The old contract declined anything outside a hand-curated root tuple, on
        the reasoning that reporting builtins "would flood the gate". True of
        reporting them DEAD — but they are ALIVE, and deciding them correctly
        floods nothing. `[M]` 2026-08-10 the curated list omitted ``mpmath``
        (the dependency the whole Peierls reference family rests on),
        ``functools``, ``dataclasses``, ``typing`` and every builtin: **86 of
        the 90** roles the tool had been dismissing as foreign.

        The lesson is general enough to be the reason this leg is parametrised
        rather than merged into one assertion: **a curated list approximating a
        decidable predicate does not announce what it is missing.**
        """
        verdict = judge(target, ("orpheus.numerics.eigenvalue",))
        assert verdict.outcome is Outcome.ALIVE, verdict

    def test_an_unknown_bare_name_is_DECLINED(self) -> None:
        """The surviving credibility leg, and the tool's honest limit.

        A bare name that is neither local, nor a builtin, nor an importable
        module cannot be told apart from "a type the reader is expected to know"
        — so it is declined rather than reported. This is the same limit that
        leaves the ``.rst`` relative bucket ungated.
        """
        verdict = judge("NoSuchTypeAnywhereInTheTree", ("orpheus.numerics.eigenvalue",))
        assert verdict.outcome is Outcome.DECLINED, verdict

    def test_a_bare_module_name_under_a_class_role_is_DECLINED(self) -> None:
        """The guard the curated list never needed.

        Computing decidability means stdlib module names become "absolute" — so
        ``:class:`array``` would resolve as the stdlib ``array`` MODULE and read
        ALIVE, for a role Sphinx leaves broken. A bare module name is legitimate
        only under ``:mod:``, so only that role takes the absolute path.
        """
        assert judge("array", (), role="class").outcome is Outcome.DECLINED
        assert judge("array", (), role="mod").outcome is Outcome.ALIVE

    def test_a_relative_target_with_no_namespace_is_DECLINED(self) -> None:
        """An ``.rst`` page has no module context, so nothing is decidable.

        This is a real residual hole, gated as a known one rather than papered
        over: ``peierls.rst`` cites ``peierls_geometry.build_volume_kernel``
        for a module retired in ``bda76faf`` and no import-based check without
        a page-level module context can see it.
        """
        assert judge("peierls_geometry.build_volume_kernel", ()).outcome is (
            Outcome.DECLINED
        )


class TestTheEmptyNamespacePackageIsNotAModule:
    """The false-ALIVE family — the one place importing is LESS informative.

    ``importlib`` materialises a PEP 420 namespace package from a bare
    directory, so ``import orpheus.<any>.<path>`` succeeds for any directory
    under ``orpheus/``. Before 2026-08-10 that made every such reference read
    alive, and the knowledge graph was right about all five shipped instances
    while this checker was wrong about all five.
    """

    RESERVED = (
        "orpheus.derivations.continuous.spectral_resolvent",
        "orpheus.derivations.continuous.pn_method",
        "orpheus.derivations.continuous.escape_probability",
        "orpheus.derivations.continuous.spectral_collocation",
        "orpheus.derivations.continuous.spn_method",
    )

    @pytest.mark.parametrize("dotted", RESERVED)
    def test_a_readme_only_directory_does_not_resolve(self, dotted: str) -> None:
        """It imports clean, and it still is not a module."""
        import importlib

        module = importlib.import_module(dotted)
        assert module.__file__ is None, (
            f"{dotted} now has a __file__ — if the reserved solver was "
            f"implemented, drop it from RESERVED and restore the live :mod: role."
        )
        assert _is_empty_namespace_package(module)
        assert resolve(dotted) == (False, "empty-namespace-package")

    def test_a_namespace_package_WITH_content_stays_alive(self) -> None:
        """The discriminator is content, not the absence of ``__init__.py``.

        ``tests/sn/`` is a namespace package carrying hundreds of modules. A
        rule keyed on ``__file__ is None`` alone would call it dead and redden
        every xref into the test tree.
        """
        import tests.sn

        assert tests.sn.__file__ is None, "premise: tests.sn is a namespace package"
        assert not _is_empty_namespace_package(tests.sn)
        assert resolve("tests.sn") == (True, None)

    def test_an_ordinary_module_is_untouched(self) -> None:
        import orpheus.numerics.face_layout as module

        assert module.__file__ is not None
        assert not _is_empty_namespace_package(module)


class TestTheTwoHalvesOfTheInterlock:
    """``#:`` scanning and ``self.X: T`` resolution are one change, not two.

    ``orpheus/numerics/face_layout.py:89`` cites
    ``:attr:`SNMesh.bc <…augmented_mesh.SNMesh.bc>``` from inside a ``#:``
    attribute-comment block, and ``SNMesh`` sets ``self.bc: dict[...] = ...``
    in ``__init__`` — which PEP 526 records **nowhere** at runtime. So widening
    the scanned surface without teaching the resolver about ``self`` annotations
    makes this citation the gate's first output, as a FALSE red whose natural
    "fix" is deleting a correct cross-reference.
    """

    CITATION = "orpheus.sn.mesh.augmented_mesh.SNMesh.bc"

    def test_the_attribute_comment_block_is_scanned_at_all(self) -> None:
        """Half one: ``ast`` discards comments, so this needs ``tokenize``."""
        path = REPO_ROOT / "orpheus" / "numerics" / "face_layout.py"
        carrying = [b for b in iter_text_blocks(path) if self.CITATION in b.text]
        assert carrying, (
            "the `#:` block citing SNMesh.bc is not being scanned — "
            "168 roles of attribute-comment prose have gone un-gated again."
        )
        assert carrying[0].namespaces == ("orpheus.numerics.face_layout",)

    def test_a_self_annotation_resolves(self) -> None:
        """Half two, and the reason it cannot be left to ``getattr``."""
        from orpheus.sn.mesh.augmented_mesh import SNMesh

        assert not hasattr(SNMesh, "bc"), (
            "premise: `bc` is per-instance. If it became a class attribute, "
            "this interlock no longer exists and this class can retire."
        )
        assert "bc" in _self_attributes(SNMesh)
        assert resolve(self.CITATION) == (True, None)

    def test_disabling_half_two_MANUFACTURES_the_false_red(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """The mutation that proves the interlock is real, not argued.

        Blind the ``self``-annotation lookup and the citation flips to DEAD.
        That is the exact wall of false reds a surface-first ordering would
        have shipped, measured rather than asserted.

        Monkeypatched in-process on purpose — never via ``git checkout``, which
        would discard uncommitted work (lessons L28).
        """
        from tools import check_docstring_xrefs as checker

        monkeypatch.setattr(checker, "_self_attributes", lambda _: frozenset())
        assert checker.resolve(self.CITATION) == (False, "missing"), (
            "blinding _self_attributes did NOT redden the citation — the "
            "interlock this class documents is not where it claims to be."
        )

    def test_a_name_never_annotated_on_self_is_not_invented(self) -> None:
        """The negative leg: the AST scan must not report everything present."""
        from orpheus.sn.mesh.augmented_mesh import SNMesh

        assert "no_such_attribute_anywhere" not in _self_attributes(SNMesh)

    def test_an_UNANNOTATED_instance_attribute_also_resolves(self) -> None:
        """The second shape, and it was a latent false red until 2026-08-10.

        ``SNMesh.mesh`` is assigned by ``MaterialMesh._init_data`` with **no**
        annotation, so neither ``getattr`` on the class nor a ``self.x: T`` scan
        can see it — yet it reads as a live ``Mesh2D`` on any instance. It was
        latent rather than active only because every citation of it happened to
        be unqualified; the first contributor to qualify one would have been
        handed a false red on a correct reference.

        Found during #346 W1 when a sub-agent's instance-level check contradicted
        my class-level probe. The class-level probe was the wrong instrument —
        the same mistake, in the same campaign, for the third time.
        """
        from orpheus.sn.mesh.augmented_mesh import SNMesh

        assert not hasattr(SNMesh, "mesh"), "premise: `mesh` is per-instance"
        assert "mesh" in _self_attributes(SNMesh)
        assert resolve("orpheus.sn.mesh.augmented_mesh.SNMesh.mesh") == (True, None)


class TestTheTargetIsReadTheWaySphinxReadsIt:
    """``extract_target`` normalisation — positive and negative legs."""

    @pytest.mark.parametrize(
        ("raw", "expected"),
        [
            ("orpheus.numerics.operator.LinearOperator", "orpheus.numerics.operator.LinearOperator"),
            ("~orpheus.numerics.operator.LinearOperator", "orpheus.numerics.operator.LinearOperator"),
            ("Display Text <orpheus.numerics.operator.LinearOperator>", "orpheus.numerics.operator.LinearOperator"),
            (".relative.path", "relative.path"),
            # The 79-column wrap: RST joins a paragraph before parsing inline
            # markup, so this is ONE target that happens to contain a newline.
            ("orpheus.numerics.trajectory.\n    ChordOracle", "orpheus.numerics.trajectory.ChordOracle"),
            ("CylinderChordOracle.\n        apply_operator", "CylinderChordOracle.apply_operator"),
        ],
    )
    def test_a_target_spelling_normalises(self, raw: str, expected: str) -> None:
        assert extract_target(raw) == expected

    @pytest.mark.parametrize(
        "raw",
        [
            "solve_galerkin_spectral_slab(c=..., d=...)",  # a signature, not a path
            "RegionMesh(n_cells, 'equal-volume')",
            "FaceLayout[tuple[int, int, str]]",  # a subscripted generic
            "",
            # ⭐ The case that corrected the implementation. A phrase separated
            # by SPACES is not a wrap, and collapsing it yields
            # ``notadottedpathatall`` — a valid identifier the whitelist happily
            # admits. Only a line BREAK is a wrap; this is the leg that says so.
            "not a dotted path at all",
            "Angular Flux",  # would collapse onto a real class name
        ],
    )
    def test_a_non_path_is_refused(self, raw: str) -> None:
        """The whitelist's teeth — and the limit of what it can do.

        Checking POSITIVELY that the result is a dotted identifier path (rather
        than blacklisting ``(``/``[``/space) is what admits a rejoined wrap
        while refusing a signature fragment or a subscripted generic. It does
        NOT by itself make rejoining safe: that is the job of joining only
        newlines. Both halves are needed and this leg pins the second one.
        """
        assert extract_target(raw) is None


class TestProseKnowsWhichNamespaceItWasWrittenIn:
    """The context derivation: module from the path, class from the AST."""

    @pytest.mark.parametrize(
        ("relative", "expected"),
        [
            ("orpheus/numerics/eigenvalue.py", "orpheus.numerics.eigenvalue"),
            ("orpheus/numerics/__init__.py", "orpheus.numerics"),
            ("tests/test_docstring_xrefs.py", "tests.test_docstring_xrefs"),
        ],
    )
    def test_the_module_comes_from_the_path(self, relative: str, expected: str) -> None:
        assert module_name_of(REPO_ROOT / relative) == expected

    def test_a_method_docstring_resolves_against_its_CLASS(self) -> None:
        """A function does not extend the namespace; a class does.

        A role in ``EigenvalueSolver.solve_fixed_source``'s docstring must
        resolve against ``…eigenvalue.EigenvalueSolver`` and then
        ``…eigenvalue`` — not against a path containing the method name.
        """
        path = REPO_ROOT / "orpheus" / "numerics" / "eigenvalue.py"
        inner = {
            block.namespaces
            for block in iter_text_blocks(path)
            if block.namespaces and len(block.namespaces) == 2
        }
        assert ("orpheus.numerics.eigenvalue.EigenvalueSolver",
                "orpheus.numerics.eigenvalue") in inner
        assert not any(
            "solve_fixed_source" in ns for namespaces in inner for ns in namespaces
        ), "a method name leaked into a namespace — roles would resolve nowhere"

    def test_an_rst_page_carries_no_namespace(self) -> None:
        """Measured premise: this project has ZERO ``currentmodule`` directives.

        So there is nothing to derive a page's module context from, and every
        relative role on a page is honestly undecidable.
        """
        page = REPO_ROOT / "docs" / "theory" / "references" / "peierls.rst"
        blocks = list(iter_text_blocks(page))
        assert len(blocks) == 1 and blocks[0].namespaces == ()
