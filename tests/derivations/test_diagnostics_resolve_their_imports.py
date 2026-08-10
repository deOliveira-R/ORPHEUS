r"""Every script under ``derivations/`` must still RESOLVE its first-party imports.

Why this gate exists (#347)
===========================

`[M]` 2026-08-09: **15 of 38** scripts in ``derivations/diagnostics/`` had
been unrunnable since the 2026-04 package restructuring, and nothing in the
tree noticed for four months. The reason is structural, not an oversight —
``derivations/`` is invisible to **every** instrument the project owns:

* ``pyproject.toml`` ships ``include = ["orpheus*"]``, so it is not packaged;
* it has no ``__init__.py`` and pytest does not collect it, so no test
  imports it;
* Sphinx renders no ``automodule`` for it, so the ``-W`` / ``-n`` build and
  the docstring-xref checker cannot see it either.

That is the same shape as the line-wrapped-role class in #346: repairing
today's instances guarantees a fresh crop at the next rename. This gate
converts a silent rot into a red test **at the rename, for the author who
made it** — the failure becomes unrepresentable as *silence*.

Not everything in there is debris. Three of the survivors are the
instrument behind something committed:

* ``diag_175_sweep_snapshot_regen.py`` regenerates ``tests/sn/sweep_ref_2g.npy``
  and carries the structurally-independent plain-Python 2-D DD oracle that
  licenses each re-baseline;
* ``diag_s69_scanmarch_vs_window_bench.py`` produced the table published in
  ``docs/theory/methods/sn/loss_representation.rst``;
* ``diag_cin_aware_split_basis_keff.py`` is imported by the subprocess worker
  in ``tests/cp/test_peierls_rank_n_protocol.py``.

What this gate CANNOT see
=========================

Resolution is static: an ``ast`` walk plus :func:`importlib.util.find_spec`.
Three import mechanisms are invisible to it, and their silence must never be
read as coverage (`vv-principles` #23 — name the readers your instrument
cannot reach):

1. **Dynamic module loading.** ``importlib.util.spec_from_file_location(...)``
   + ``exec_module``. This was not hypothetical: ``diag_phase_g_step2_mesh_scaling.py``
   resolved every ``import`` statement it had while being broken through a
   sibling it ``exec``'d. Both were retired at #347, so the tree carries
   **zero** instances today.
2. **Imports inside a subprocess-worker string.** A ``textwrap.dedent`` body
   handed to ``python -c`` is plain text. The one live instance is in
   ``tests/cp/test_peierls_rank_n_protocol.py``, where the consuming test
   *executes* it — so that path is gated by its own test, not by this one.
   (This class already cost ten weeks once: ``15486f66`` mass-deleted
   ``diag_cin_aware_split_basis_keff`` while that worker consumed it.)
3. **A computed module name** — ``importlib.import_module(f"...{x}")``.

And resolution is not execution: a script whose imports resolve can still
fail on a changed *signature* (``TimedFullField.zeros(bulk=...)`` →
``interior=``) or a renamed attribute reached through an object rather than
an import. That second layer is genuinely out of reach for a ~1 s gate; the
scripts that must actually RUN are the three named above, and each is either
run by hand before it is relied on or exercised by its consuming test.

Decoder discipline
==================

The flag fires on three distinct states, each of which needs its own
control rather than one control for the instrument
(`vv-principles` #8, the ninth-class METHOD WARNING — a predicate written
to make one decision carries its other meanings along):

===================================  ================================================
state                                control in :class:`TestTheGateItselfHasTeeth`
===================================  ================================================
module does not exist                ``test_it_flags_a_dead_module``
module exists, import raises         ``test_it_flags_a_module_that_raises_on_import``
module exists, attribute missing     ``test_it_flags_a_missing_attribute``
===================================  ================================================

The known false-positive source is a module-level ``__getattr__`` (PEP 562),
which serves names that ``hasattr`` on a fresh import may not show; it is
excluded explicitly and pinned by ``test_a_module_level_getattr_is_not_flagged``.
"""

from __future__ import annotations

import ast
import importlib
import importlib.util
import pathlib

import pytest

pytestmark = pytest.mark.foundation

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
DERIVATIONS_ROOT = REPO_ROOT / "derivations"

#: Import roots this gate resolves. Third-party imports are the venv's problem
#: (a missing one fails loudly at install time); these two are the ones a
#: rename inside THIS repository can silently orphan.
FIRST_PARTY = ("orpheus", "tests")


def _scripts() -> list[pathlib.Path]:
    return sorted(
        p
        for p in DERIVATIONS_ROOT.rglob("*.py")
        if "__pycache__" not in p.parts
    )


def _is_first_party(dotted: str) -> bool:
    return dotted.split(".")[0] in FIRST_PARTY


def _unresolvable(source: str, *, where: str) -> list[str]:
    """Every first-party name in ``source`` that no longer resolves.

    Returns one human-readable line per defect, empty when the file is clean.
    """
    problems: list[str] = []
    for node in ast.walk(ast.parse(source, filename=where)):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if _is_first_party(alias.name) and not _resolves(alias.name):
                    problems.append(f"line {node.lineno}: no module {alias.name!r}")
        elif isinstance(node, ast.ImportFrom):
            module = node.module or ""
            if node.level or not _is_first_party(module):
                continue  # relative import — derivations/ is not a package
            if not _resolves(module):
                problems.append(f"line {node.lineno}: no module {module!r}")
                continue
            try:
                imported = importlib.import_module(module)
            except Exception as exc:  # noqa: BLE001 — any import-time failure counts
                problems.append(
                    f"line {node.lineno}: importing {module!r} raised "
                    f"{type(exc).__name__}: {exc}"
                )
                continue
            if hasattr(imported, "__getattr__"):
                continue  # PEP 562 lazy surface — hasattr is not decisive
            for alias in node.names:
                if alias.name == "*":
                    continue
                if hasattr(imported, alias.name):
                    continue
                if _resolves(f"{module}.{alias.name}"):
                    continue  # a submodule, not an attribute
                problems.append(
                    f"line {node.lineno}: {module!r} has no {alias.name!r}"
                )
    return problems


def _resolves(dotted: str) -> bool:
    try:
        return importlib.util.find_spec(dotted) is not None
    except (ImportError, ModuleNotFoundError, ValueError):
        return False


@pytest.mark.parametrize(
    "script", _scripts(), ids=lambda p: str(p.relative_to(DERIVATIONS_ROOT))
)
def test_every_derivations_script_resolves_its_first_party_imports(
    script: pathlib.Path,
) -> None:
    """A rename that orphans a script reds HERE, not four months later."""
    problems = _unresolvable(script.read_text(), where=str(script))
    assert not problems, (
        f"{script.relative_to(REPO_ROOT)} no longer resolves:\n  "
        + "\n  ".join(problems)
        + "\n\nA rename orphaned it. Decide per #347's discriminator — does "
        "anything OUTSIDE derivations/ depend on this file?\n"
        "  * a test / doc / production module names it as the instrument "
        "behind a pinned baseline, a published number, or a pointer meant to "
        "be followed  -> REPOINT it;\n"
        "  * it calls itself a gate and no committed test covers its claim  "
        "-> PROMOTE it into tests/;\n"
        "  * nothing depends on it and its finding has landed  -> RETIRE it "
        "(git keeps the archaeology).\n"
        "Do NOT blanket-repoint: a script that imports cleanly while still "
        "computing against retired APIs is strictly worse than an honest "
        "ImportError."
    )


class TestTheGateItselfHasTeeth:
    """One control per state the predicate maps to True, not one per instrument.

    Without these the gate is exactly the shape `vv-principles` #17 warns
    about — an all-green verdict from an instrument nobody proved can fire.
    """

    def test_it_flags_a_dead_module(self) -> None:
        problems = _unresolvable(
            "from orpheus.sn.quadrature import Quadrature\n", where="<probe>"
        )
        assert problems == ["line 1: no module 'orpheus.sn.quadrature'"], problems

    def test_it_flags_a_missing_attribute(self) -> None:
        problems = _unresolvable(
            "from orpheus.sn.loss_representation import transport_sweep\n",
            where="<probe>",
        )
        assert len(problems) == 1 and "transport_sweep" in problems[0], problems

    def test_it_flags_a_module_that_raises_on_import(self, monkeypatch) -> None:
        """The middle state — resolvable, but dead on arrival.

        Distinct from a missing module: ``find_spec`` says yes and the
        failure only appears at execution, which is precisely the state a
        ``find_spec``-only scan reports as clean.
        """
        boom = RuntimeError("deliberate")

        def explode(name: str, *args, **kwargs):
            if name == "orpheus.numerics.convergence":
                raise boom
            return importlib.__dict__["__import_module_real"](name, *args, **kwargs)

        monkeypatch.setattr(
            importlib, "__import_module_real", importlib.import_module, raising=False
        )
        monkeypatch.setattr(importlib, "import_module", explode)

        problems = _unresolvable(
            "from orpheus.numerics.convergence import IterationRecord\n",
            where="<probe>",
        )
        assert len(problems) == 1, problems
        assert "raised RuntimeError: deliberate" in problems[0], problems

    def test_a_module_level_getattr_is_not_flagged(self) -> None:
        """The known false positive, pinned so a future tightening cannot

        silently re-introduce it. ``typing`` defines a module-level
        ``__getattr__``, so an absent name there must NOT be reported.
        """
        assert hasattr(importlib.import_module("typing"), "__getattr__")
        assert (
            _unresolvable(
                "from typing import a_name_typing_does_not_export\n", where="<probe>"
            )
            == []
        )

    def test_a_live_import_is_not_flagged(self) -> None:
        """The negative leg: the gate must not cry wolf on healthy imports."""
        assert (
            _unresolvable(
                "import orpheus.sn.solver\n"
                "from orpheus.numerics.quadrature import Quadrature\n"
                "from tests.sn._test_helpers import sweep_once\n",
                where="<probe>",
            )
            == []
        )

    def test_the_gate_actually_covers_the_tree(self) -> None:
        """A parametrize over an empty glob is green and gates nothing.

        The ``rglob`` runs at collection time against a path built from
        ``__file__``; if ``derivations/`` moves, the row count silently
        drops to zero rather than reddening (`vv-principles` #8's
        signature-tautological class).
        """
        scripts = _scripts()
        assert len(scripts) >= 20, (
            f"only {len(scripts)} scripts found under {DERIVATIONS_ROOT} — "
            "either the tree moved (repoint DERIVATIONS_ROOT) or a mass "
            "retirement happened; a near-empty parametrize gates nothing."
        )
