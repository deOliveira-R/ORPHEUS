r"""The new condition family shares no qualname with the retiring tree (wave 2, WF2.3 slice S3, gate 8 of 6.3).

Spec: ``.claude/plans/wave2_boundary_carve.md`` section 2.7 (why the diffusion type is ``ZeroFluxCondition`` and not
``ZeroFluxBoundary``), decisions 29 and 30, and the finding of section 0: ``tests/gates/values/test_content_identity_relocation.py``'s
``test_no_two_nameable_orpheus_classes_share_a_qualname`` fails when two content-keyed classes share a qualname, because a
content key names an orpheus type by its qualname. The ruling of the thirteenth exchange is that the gate takes NO exemption (an
exemption would be a dual spelling under another name), so the carve moves classes in under names that are fresh. The old
``geometry/boundary/zero_flux.py`` class ``ZeroFluxBoundary`` cannot die before WF2.4 (it is named by 8 production and 10 test
files of slices S4 and S5), so the diffusion type minted in S3 is named ``ZeroFluxCondition`` and lives in
``orpheus/diffusion/zero_flux.py``.

The claim, by an AST pass (never a grep: ``instrument-doctrine`` X2): over ``orpheus/``,

* the **new family** is every class defined in the conditions package ``orpheus/transport/boundary_conditions/`` and in
  ``orpheus/diffusion/zero_flux.py``;
* the **retiring tree** is every class defined in ``orpheus/geometry/boundary/`` and in ``orpheus/transport/boundary_realization.py``
  (the two homes WF2.4 deletes);
* no qualname (the dotted class path through nesting) is in both, EXCEPT exception classes. The exclusion is stated, because the
  plan states it (section 0, decision 29): ``BoundaryError``, ``SubmarkovViolationError`` and ``BoundaryResponseNotPositiveError``
  are redefined by ``errors.py`` and die in the old tree with their package at WF2.4; the qualname gate covers content-keyed
  classes (a ``ContentIdentity`` subclass, a frozen dataclass, an ``Enum``), not exceptions. An exception is a class one of whose
  bases is a builtin exception or an exception of the scanned set.

Why red at the base tip (``0757a8d9``): ``orpheus/diffusion/zero_flux.py`` does not exist, so the new family has no
``ZeroFluxCondition`` and the gate's own population check (the type is minted in S3 and is the only class that the S3 file defines
by the plan) fails. The controls and the retiring tree's classes expire with WF2.4, which deletes the old tree: the gate is retired
there with the classes it compares against. The scanner rows (``test_the_scanner_sees_what_it_must``) are green today and are the positive controls of the
red row: they run on planted trees under ``tmp_path`` and need no production module.

Mutations for the reviewers: the type named ``ZeroFluxBoundary`` (a copy of the new module with the rename: the planted control
does exactly this and reads one shared qualname); a class ``AlbedoBoundary`` or ``WhiteBoundary`` defined in the conditions
package; a nested class that repeats an old nested qualname.
"""

from __future__ import annotations

import ast
import builtins
from pathlib import Path

import pytest

import orpheus as _orpheus

pytestmark = pytest.mark.foundation

_THIS = "tests/gates/diffusion/test_condition_family_qualnames.py"
ROOT = Path(_orpheus.__file__).parent.parent

#: The homes of the new family and of the retiring tree, relative to the tree root.
NEW_PACKAGE = "orpheus/transport/boundary_conditions"
NEW_FILES = ("orpheus/diffusion/zero_flux.py",)
OLD_PACKAGE = "orpheus/geometry/boundary"
OLD_FILES = ("orpheus/transport/boundary_realization.py",)


def class_table(path: Path) -> dict[str, list[str]]:
    """``{qualname: [base names]}`` of every class defined in ``path`` (an AST walk; a class in a function body reads
    ``f.<locals>.C``, as ``__qualname__`` does)."""
    found: dict[str, list[str]] = {}

    def walk(node: ast.AST, stack: tuple[str, ...]) -> None:
        for child in ast.iter_child_nodes(node):
            if isinstance(child, ast.ClassDef):
                found[".".join((*stack, child.name))] = [ast.unparse(base) for base in child.bases]
                walk(child, (*stack, child.name))
            elif isinstance(child, ast.FunctionDef | ast.AsyncFunctionDef):
                walk(child, (*stack, child.name, "<locals>"))
            else:
                walk(child, stack)

    walk(ast.parse(path.read_text()), ())
    return found


def gather(root: Path, package: str, files: tuple[str, ...]) -> dict[str, list[str]]:
    """The classes of every ``.py`` under ``root/package`` and of each of ``files``, merged by qualname (the bases of a homonym
    are concatenated). A listed file that does not exist is skipped here and asserted present by the gate that needs it (the
    controls must run at a base tip where the S3 file is not yet minted)."""
    paths = sorted((root / package).rglob("*.py")) + [root / f for f in files if (root / f).exists()]
    merged: dict[str, list[str]] = {}
    for path in paths:
        for qualname, bases in class_table(path).items():
            merged.setdefault(qualname, []).extend(bases)
    return merged


def exceptions_of(table: dict[str, list[str]]) -> set[str]:
    """The qualnames of ``table`` whose bases reach a builtin exception, transitively within ``table``."""
    flagged: set[str] = set()
    changed = True
    while changed:
        changed = False
        for qualname, bases in table.items():
            if qualname in flagged:
                continue
            for base in bases:
                builtin = getattr(builtins, base, None)
                if (isinstance(builtin, type) and issubclass(builtin, BaseException)) or base in flagged:
                    flagged.add(qualname)
                    changed = True
                    break
    return flagged


def shared_qualnames(root: Path) -> tuple[set[str], set[str], int, int]:
    """``(shared, shared_with_exceptions, n_new, n_old)``: the qualnames in both families, with and without the exception exclusion."""
    new = gather(root, NEW_PACKAGE, NEW_FILES)
    old = gather(root, OLD_PACKAGE, OLD_FILES)
    both = set(new) & set(old)
    excluded = exceptions_of(new) | exceptions_of(old)
    return both - excluded, both, len(new), len(old)


def test_the_scanner_sees_what_it_must(tmp_path: Path) -> None:
    r"""The positive controls, on planted trees (green today; reddens on a scanner that reads nothing or reports everything).

    * the REAL tree: the retiring tree holds the known members ``ZeroFluxBoundary``, ``AlbedoBoundary`` and ``WhiteBoundary``
      (the filter finds a member of each file shape), and the unfiltered overlap holds exactly the three exceptions the plan names
      (so the exclusion is exercised: it removes something real);
    * a planted new module that defines ``ZeroFluxBoundary`` beside a planted old one reads ONE shared qualname, while the same
      module with the class named ``ZeroFluxCondition`` reads none;
    * a nested class repeats a nested qualname only when its whole path repeats (``A.Inner`` vs ``B.Inner`` is not shared);
    * a planted shared exception name is excluded and a planted shared dataclass name is not.
    """
    old = gather(ROOT, OLD_PACKAGE, OLD_FILES)
    assert {"ZeroFluxBoundary", "AlbedoBoundary", "WhiteBoundary"} <= set(old), sorted(old)[:8]
    shared, both, n_new, n_old = shared_qualnames(ROOT)
    assert both >= {"BoundaryError", "SubmarkovViolationError", "BoundaryResponseNotPositiveError"}
    assert n_new > 20 and n_old > 20, (n_new, n_old)

    def plant(new_text: str, old_text: str) -> Path:
        for rel, text in ((f"{NEW_PACKAGE}/m.py", new_text), (NEW_FILES[0], ""), (f"{OLD_PACKAGE}/m.py", old_text), (OLD_FILES[0], "")):
            target = tmp_path / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(text)
        return tmp_path

    root = plant("class ZeroFluxBoundary:\n    pass\n", "class ZeroFluxBoundary:\n    pass\n")
    assert shared_qualnames(root)[0] == {"ZeroFluxBoundary"}
    root = plant("class ZeroFluxCondition:\n    pass\n", "class ZeroFluxBoundary:\n    pass\n")
    assert shared_qualnames(root)[0] == set()
    root = plant("class A:\n    class Inner:\n        pass\n", "class B:\n    class Inner:\n        pass\n")
    assert shared_qualnames(root)[0] == set()
    root = plant("class A:\n    class Inner:\n        pass\n", "class A:\n    class Inner:\n        pass\n")
    assert shared_qualnames(root)[0] == {"A", "A.Inner"}
    root = plant("class E(ValueError):\n    pass\n\nclass D:\n    pass\n", "class E(ValueError):\n    pass\n\nclass D:\n    pass\n")
    shared, both, _, _ = shared_qualnames(root)
    assert both == {"E", "D"} and shared == {"D"}


@pytest.mark.rests_on(f"{_THIS}::test_the_scanner_sees_what_it_must")
def test_the_condition_family_shares_no_qualname_with_the_old_tree() -> None:
    r"""Gate 8 of 6.3. The new family (the conditions package and ``orpheus/diffusion/zero_flux.py``) defines ``ZeroFluxCondition``
    and shares no non-exception qualname with the retiring tree (``orpheus/geometry/boundary/`` and
    ``orpheus/transport/boundary_realization.py``). The population is printed in the failure: classes of each side.

    Red input at the base tip: ``orpheus/diffusion/zero_flux.py`` does not exist (the carve mints it). Mutation: the type named
    ``ZeroFluxBoundary``; any class of the old tree's names defined in the new family.
    """
    zero_flux = ROOT / NEW_FILES[0]
    assert zero_flux.exists(), f"{NEW_FILES[0]} does not exist: slice S3 mints ZeroFluxCondition there"
    assert "ZeroFluxCondition" in class_table(zero_flux), f"{NEW_FILES[0]} does not define ZeroFluxCondition: {sorted(class_table(zero_flux))}"
    shared, both, n_new, n_old = shared_qualnames(ROOT)
    assert not shared, f"{len(shared)} qualnames shared between the new family ({n_new} classes) and the retiring tree ({n_old}): {sorted(shared)}"
    assert "ZeroFluxBoundary" in gather(ROOT, OLD_PACKAGE, OLD_FILES) and "ZeroFluxBoundary" not in class_table(zero_flux)
