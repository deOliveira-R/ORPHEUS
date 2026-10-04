r"""Step 1 of #405 P3: the manifest — what ran, pinned by normalised AST, distribution and interpreter.

Spec ``.claude/plans/reference_p3_spec.md`` §1.1, gates M1.1–M1.7. Every production name is resolved through
``_traced_memo_api`` at run time; every mutated source is a copy in ``tmp_path`` (``_traced_memo_synthetic``).
"""
from __future__ import annotations

import json
import os
import site
import sysconfig
import textwrap

import numpy as np
import pytest

from . import _traced_memo_api as api
from ._traced_memo_synthetic import Package

pytestmark = pytest.mark.foundation


@pytest.fixture
def package(tmp_path):
    pkg = Package(tmp_path)
    with api.cache_root(tmp_path / "cache"):
        yield pkg
    pkg.close()


SOURCE = textwrap.dedent('''
    """Module doc."""
    import math

    LIMIT = 3


    def target(x, y=2):
        """Doc one."""
        # a comment
        total = x + y
        return total * LIMIT


    def neighbour(x):
        return x - 1


    class Holder:
        WIDTH = 4

        def act(self, x):
            return x + self.WIDTH
''')


def _edited(old: str, new: str) -> str:
    if SOURCE.count(old) != 1:
        raise AssertionError(f"edit target {old!r} occurs {SOURCE.count(old)} times")
    return SOURCE.replace(old, new)


#: (edit name, old, new, the def's digest moves, the skeleton moves)
EDITS = [
    ("docstring", '"""Doc one."""', '"""A different doc."""', False, False),
    ("comment", "# a comment", "# another comment entirely", False, False),
    ("blank-lines-above", "import math\n", "import math\n\n\n\n", False, False),
    ("module-docstring", '"""Module doc."""', '"""Another module doc."""', False, False),
    ("reformat", "total = x + y", "total = (x\n                 + y)", False, False),
    ("body-operator", "total = x + y", "total = x - y", True, False),
    ("body-constant", "return total * LIMIT", "return total * LIMIT * 1", True, False),
    ("signature-default", "def target(x, y=2):", "def target(x, y=3):", True, True),
    ("decorator", "def target(x, y=2):", "@staticmethod\ndef target(x, y=2):", True, True),
    ("neighbour-body", "return x - 1", "return x - 2", False, False),
    ("module-constant", "LIMIT = 3", "LIMIT = 4", False, True),
    ("import", "import math", "import math\nimport cmath", False, True),
    ("def-added", "class Holder:", "def added():\n    return 0\n\n\nclass Holder:", False, True),
    ("class-attribute", "WIDTH = 4", "WIDTH = 5", False, True),
    ("method-body", "return x + self.WIDTH", "return x + self.WIDTH + 0", False, False),
]


@pytest.mark.parametrize("edit,old,new,def_moves,skeleton_moves", EDITS, ids=[e[0] for e in EDITS])
def test_m1_1_and_m1_2_the_digests_move_exactly_with_what_can_change_an_answer(edit, old, new, def_moves, skeleton_moves):
    """M1.1 (the def's digest) and M1.2 (the module skeleton): docstrings, comments and layout never move a
    digest; a body edit moves its def's and nobody else's; a top-level statement, a signature, a decorator, an
    import or a def added moves the skeleton. Each row is one edit and asserts BOTH digests, so a skeleton that
    kept function bodies, or a def digest that kept docstrings, reds a named row."""
    after = _edited(old, new)
    before_def, after_def = api.function_digest(SOURCE, "target"), api.function_digest(after, "target")
    if before_def is None or after_def is None:
        raise AssertionError("function_digest did not find 'target'")
    assert (before_def != after_def) is def_moves, f"{edit}: the def digest moved={before_def != after_def}"
    before_skel, after_skel = api.skeleton_digest(SOURCE), api.skeleton_digest(after)
    assert (before_skel != after_skel) is skeleton_moves, f"{edit}: the skeleton moved={before_skel != after_skel}"


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_manifest.py::test_m1_1_and_m1_2_the_digests_move_exactly_with_what_can_change_an_answer')
def test_m1_3a_a_method_and_the_last_def_of_a_name_resolve_by_qualname():
    """M1.3 (a): ``Holder.act`` is found by qualname, and a second ``def target`` shadows the first, as Python
    binds it: the digest is the LAST def's, so a validation reads the code that would run."""
    shadowed = SOURCE + "\n\ndef target(x, y=2):\n    return 0\n"
    alone = "def target(x, y=2):\n    return 0\n"
    assert api.function_digest(SOURCE, "Holder.act") is not None
    assert api.function_digest(SOURCE, "Holder.nope") is None
    assert api.function_digest(shadowed, "target") == api.function_digest(alone, "target")


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_manifest.py::test_m1_3a_a_method_and_the_last_def_of_a_name_resolve_by_qualname')
def test_m1_3b_nested_code_maps_to_its_outermost_def(package):
    """M1.3 (b): a nested def, a lambda and a generator expression inside ``nested_user`` are pinned by
    ``nested_user``'s digest; a method by ``Box.method``; no ``<locals>`` qualname reaches the manifest."""
    alpha = package.module("alpha")
    _, manifest = api.trace_call(alpha.generate.__wrapped__, 1.5)
    names = {q for rel, q, *_ in manifest.functions if rel == package.relpath("alpha")}
    assert {"nested_user", "Box.method", "helper", "generate", "_count"} <= names, names
    assert not any("<" in q for _, q, *_ in manifest.functions), manifest.functions


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_manifest.py::test_m1_3b_nested_code_maps_to_its_outermost_def')
def test_m1_3c_signature_scopes_are_pinned_by_the_skeleton_not_by_a_body(package):
    """M1.3 (c), Python 3.14: evaluating a def's lazy annotations runs a generated ``__annotate__`` code object
    reported at the def's line, and a PEP 695 generic runs a ``<generic parameters of ...>`` scope. Neither is
    refused, and neither pins the def's BODY: ``annotated`` (whose annotations were read, never its body) is
    absent from the manifest while the module's skeleton, which keeps every signature, is present. Pinning the
    body here would make a parent stale through its own manifest when only its child's body changed, so the
    recursive witness M3.12 (b) would pass for the wrong reason (``[M]`` 2026-10-04: the prototype's first
    run listed the child def ``inner`` in the parent's manifest, through ``get_type_hints``)."""
    delta = package.module("delta")
    value, manifest = api.trace_call(delta.reads_annotations, 2.0)
    rel = package.relpath("delta")
    names = {q for r, q, *_ in manifest.functions if r == rel}
    assert names == {"reads_annotations", "generic", "Holder.scaled"}, names
    assert rel in {r for r, _ in manifest.modules}
    assert value == 5.0  # 2.0 + three hints (x, y, return)


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_manifest.py::test_m1_3b_nested_code_maps_to_its_outermost_def')
def test_m1_4_the_manifest_is_what_ran(package):
    """M1.4: the traced defs are exactly the ones that ran (``unused`` and ``Box.other`` absent: the negative
    leg that a whole-module hash would fail), the module's skeleton is present, numpy is pinned by its installed
    version, no stdlib file is a source entry, and the dataclass's generated ``__init__`` (``<string>``) is
    dropped while the dataclass definition is pinned by the skeleton."""
    alpha = package.module("alpha")
    _, manifest = api.trace_call(alpha.solve.__wrapped__, 4)
    _, manifest_g = api.trace_call(alpha.generate.__wrapped__, 1.5)
    rel = package.relpath("alpha")
    ran = {q for r, q, *_ in manifest_g.functions if r == rel}
    assert "unused" not in ran and "Box.other" not in ran, ran
    assert {"solve", "_count"} <= {q for r, q, *_ in manifest.functions if r == rel}
    assert rel in {r for r, _ in manifest.modules}
    assert dict(manifest_g.distributions).get("numpy") == np.__version__, manifest_g.distributions
    stdlib = os.path.realpath(sysconfig.get_paths()["stdlib"])
    assert not any(r.startswith(stdlib) or r.startswith("json/") for r, *_ in manifest.functions)
    assert not any(r.startswith("<") for r, *_ in manifest.functions)
    assert len(manifest.functions) >= 2 and len(manifest.modules) >= 1  # the population is not empty


def test_m1_5_every_shape_of_traced_file_is_classified():
    """M1.5: a numpy file is the ``numpy`` distribution, a stdlib file the interpreter, frozen and generated
    code is dropped, a source file is hashed, and a site-packages file of NO distribution is refused
    (``Unpinnable``): nothing that ran escapes the manifest unpinned."""
    assert api.classify(np.__file__) == ("distribution", "numpy")
    assert api.classify(json.__file__)[0] == "python"
    assert api.classify("<frozen importlib._bootstrap>")[0] == "dropped"
    assert api.classify("<string>")[0] == "dropped"
    assert api.classify(__file__)[0] == "source"
    orphan = os.path.join(site.getsitepackages()[0], "no_distribution_owns_this_xyz", "mod.py")
    with pytest.raises(api.name("Unpinnable"), match="no installed distribution"):
        api.classify(orphan)


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_manifest.py::test_m1_3a_a_method_and_the_last_def_of_a_name_resolve_by_qualname')
def test_m1_6_code_with_no_def_at_its_line_is_refused(package, tmp_path):
    """M1.6: a function compiled from a string under a REAL file's name, at a line where that file defines
    no such function, has no source the validation could re-hash: the manifest refuses it."""
    alpha_path = str(package.path("alpha"))
    namespace: dict = {}
    exec(compile("\n\ndef phantom(x):\n    return x * 2\n", alpha_path, "exec"), namespace)
    with pytest.raises(api.name("Unpinnable"), match="no def of that name"):
        api.trace_call(namespace["phantom"], 1.0)


def test_m1_7_the_interpreter_is_in_the_manifest(package):
    """M1.7: the manifest names the interpreter (version and cache tag), which pins the stdlib and the
    ``<string>`` code it generates; the platform tag names the OS, the machine, the optimisation level and
    the BLAS (the key, §3)."""
    import sys

    alpha = package.module("alpha")
    _, manifest = api.trace_call(alpha.helper, 1.0)
    identity = api.interpreter_identity(manifest)
    assert sys.version in identity and sys.implementation.cache_tag in identity
    tag = api.platform_tag()
    assert sys.platform in tag and f"O{sys.flags.optimize}" in tag
