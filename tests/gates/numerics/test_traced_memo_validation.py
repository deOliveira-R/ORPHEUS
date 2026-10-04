r"""Step 2 of #405 P3: validation against the checkout, and the payload codec.

Spec ``.claude/plans/reference_p3_spec.md`` §1.2, gates M2.1–M2.9. The manifests here come from ``trace_call``
in this process on a fresh synthetic package (step 1); every edit is to that package's copy in ``tmp_path``.
"""
from __future__ import annotations

import dataclasses
import importlib.metadata
import math
import os
import random
import struct
import subprocess
import sys
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


def _generate_manifest(package):
    alpha = package.module("alpha")
    value, manifest = api.trace_call(alpha.generate.__wrapped__, 1.5)
    assert api.validate(manifest) == (), "a fresh manifest validates"  # the precondition of every row
    return manifest


#: (row, module, old, new, stale?, the fragment the reason names)
WITNESSES = [
    ("traced-body", "alpha", "return x * SCALE", "return x * SCALE * 1.0", True, "helper"),
    ("traced-method", "alpha", "return x + 2.0", "return x + 2.5", True, "Box.method"),
    ("traced-nested-lambda", "alpha", "(lambda z: z * 0.5)", "(lambda z: z * 0.25)", True, "nested_user"),
    ("constant-read", "alpha", "SCALE = 3.0", "SCALE = 3.5", True, "skeleton"),
    ("constant-unread", "alpha", "UNREAD = 7.0", "UNREAD = 8.0", True, "skeleton"),
    ("untraced-body", "alpha", "return x - 1.0", "return x - 9.0", False, ""),
    ("untraced-method", "alpha", "return x - 2.0", "return x - 9.0", False, ""),
    ("docstring", "alpha", '"""Doc of a traced helper."""', '"""Changed."""', False, ""),
    ("comment-added", "alpha", "def helper(x):", "# a comment\ndef helper(x):", False, ""),
    ("module-never-imported", "gamma", "NEVER = 1.0", "NEVER = 2.0", False, ""),
]


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_manifest.py::test_m1_1_and_m1_2_the_digests_move_exactly_with_what_can_change_an_answer', 'tests/gates/numerics/test_traced_memo_manifest.py::test_m1_4_the_manifest_is_what_ran')
@pytest.mark.parametrize("row,module,old,new,stale,fragment", WITNESSES, ids=[w[0] for w in WITNESSES])
def test_m2_1_validation_misses_exactly_when_something_that_ran_changed(package, row, module, old, new, stale, fragment):
    """M2.1, the X1 witnesses at the manifest: an edit inside a traced def (a body, a method, a lambda nested in
    a def) or to a top-level statement of a module whose body ran (a constant, read or not: the skeleton is
    per module, so ``constant-unread`` is a DECLARED spurious miss, sound and priced) is a reason naming what
    changed; an edit to an untraced def of the same module, a docstring, a comment, or a module the call
    never imported, is not."""
    manifest = _generate_manifest(package)
    package.edit(module, old, new)
    reasons = api.validate(manifest)
    assert bool(reasons) is stale, f"{row}: reasons={reasons}"
    if stale:
        assert any(fragment in r for r in reasons), f"{row}: no reason names {fragment!r}: {reasons}"


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_validation.py::test_m2_1_validation_misses_exactly_when_something_that_ran_changed')
def test_m2_2_a_removed_def_and_a_removed_file_are_stale(package):
    """M2.2: a traced def that no longer exists, and a module file that is gone, are reasons (a lookup by
    qualname that returns nothing must not read as "unchanged")."""
    manifest = _generate_manifest(package)
    package.edit("alpha", "def helper(x):", "def helper_renamed(x):")
    assert any("helper" in r for r in api.validate(manifest))
    os.remove(package.path("alpha"))
    assert any("gone" in r for r in api.validate(manifest))


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_manifest.py::test_m1_7_the_interpreter_is_in_the_manifest')
def test_m2_3_a_distribution_version_and_the_interpreter_are_validated(package):
    """M2.3: a manifest recording another numpy version, or another interpreter, is stale. The rows edit the
    MANIFEST (the installed environment is not mutable from a gate); the reason names the distribution."""
    manifest = _generate_manifest(package)
    distributions = tuple((n, "0.0.0-not-installed" if n == "numpy" else v) for n, v in manifest.distributions)
    assert any("numpy" in r for r in api.validate(dataclasses.replace(manifest, distributions=distributions)))
    assert any("python" in r for r in api.validate(dataclasses.replace(manifest, python="CPython 2.7")))
    gone = manifest.distributions + (("no-such-distribution-xyz", "1.0"),)
    assert any("no-such-distribution-xyz" in r for r in api.validate(dataclasses.replace(manifest, distributions=gone)))


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_validation.py::test_m2_1_validation_misses_exactly_when_something_that_ran_changed')
def test_m2_4_validation_imports_and_runs_nothing(package, tmp_path):
    """M2.4: a manifest of the ``marker`` module (whose import writes a file) validates in a FRESH interpreter
    without the marker reappearing and without importing the package; the positive control imports it and the
    marker appears. Validation reads files; it never executes them."""
    marker_module = package.module("marker")
    _, manifest = api.trace_call(marker_module.plain.__wrapped__, 1.0)
    package.marker.unlink()
    probe = tmp_path / "probe.py"
    probe.write_text(textwrap.dedent(f'''
        import pickle, sys
        sys.path[:0] = {[str(package.root), os.getcwd()]!r}
        from {api.MODULE} import validate
        manifest = pickle.loads(open({str(tmp_path / "m.pkl")!r}, "rb").read())
        reasons = validate(manifest)
        print("REASONS", reasons)
        print("IMPORTED", any(m.startswith({package.name!r}) for m in sys.modules))
        if len(sys.argv) > 1:
            import {package.name}.marker
    '''))
    import pickle

    (tmp_path / "m.pkl").write_bytes(pickle.dumps(manifest))
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([os.getcwd(), *sys.path]))
    run = subprocess.run([sys.executable, "-O", str(probe)], capture_output=True, text=True, env=env, cwd=str(tmp_path))
    assert "REASONS ()" in run.stdout and "IMPORTED False" in run.stdout, run.stdout + run.stderr
    assert not package.marker.exists(), "validation ran the module body"
    control = subprocess.run([sys.executable, "-O", str(probe), "import"], capture_output=True, text=True, env=env, cwd=str(tmp_path))
    assert package.marker.exists(), "the control: importing the module writes the marker " + control.stderr


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_validation.py::test_m2_1_validation_misses_exactly_when_something_that_ran_changed')
def test_m2_5_a_later_edit_in_the_same_process_is_seen(package):
    """M2.5: validation has no in-process cache that outlives an edit: validate (valid), edit, validate again
    in the same process (stale), restore the text (valid again). A cache keyed by path, or by modification
    time at coarse resolution, would serve the first answer."""
    manifest = _generate_manifest(package)
    original = package.path("alpha").read_text()
    package.edit("alpha", "return x * SCALE", "return x * SCALE * 2.0")
    assert api.validate(manifest) != ()
    package.path("alpha").write_text(original)
    assert api.validate(manifest) == ()


# ── the payload codec ──────────────────────────────────────────────────────────


def _bits(x: float) -> bytes:
    return struct.pack("<d", x)


def test_m2_6_floats_cross_the_payload_bit_for_bit():
    """M2.6: every finite double round-trips bit for bit through the JSON payload: a hand population
    (``-0.0``, the smallest subnormal, the largest double, ``1 + ulp``, ``1/3``, ``0.1``) and 20 000 random bit
    patterns (seed 20261004, NaN and infinities excluded: no reading holds one). Python ``int`` beyond 2**53
    and ``bool`` keep their type and value."""
    rng = random.Random(20261004)
    population = [-0.0, 0.0, 5e-324, -5e-324, sys.float_info.max, -sys.float_info.max, math.nextafter(1.0, 2.0),
                  1.0 / 3.0, 0.1, 2.2250738585072014e-308]
    while len(population) < 20_010:
        x = struct.unpack("<d", rng.getrandbits(64).to_bytes(8, "little"))[0]
        if math.isfinite(x):
            population.append(x)
    tree, arrays = api.encode_payload(tuple(population))
    back = api.decode_payload(tree, arrays, set())
    assert len(back) == len(population)
    bad = [i for i, (a, b) in enumerate(zip(population, back)) if _bits(a) != _bits(b) or type(b) is not float]
    assert not bad, f"{len(bad)} of {len(population)} floats moved, first at {bad[:3]}"
    for value in (2**60 + 1, True, False, None, "text", 7):
        tree, arrays = api.encode_payload((value,))
        (out,) = api.decode_payload(tree, arrays, set())
        assert out == value and type(out) is type(value), (value, out)


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_validation.py::test_m2_6_floats_cross_the_payload_bit_for_bit')
def test_m2_7_arrays_keep_dtype_shape_and_bits_and_load_fresh_and_read_only():
    """M2.7: an array keeps its dtype, shape and bits; every decode hands out a FRESH array (no shared memory
    with the last load) that is READ-ONLY, so a consumer that writes into a served payload raises instead of
    corrupting a shared hit. A numpy scalar keeps its numpy type."""
    rng = np.random.default_rng(20261004)
    value = (rng.standard_normal((3, 4, 2)), np.arange(5, dtype=np.int64), np.array([True, False]), np.float64(1.25))
    tree, arrays = api.encode_payload(value)
    first = api.decode_payload(tree, arrays, set())
    second = api.decode_payload(tree, arrays, set())
    for original, a, b in zip(value[:3], first[:3], second[:3]):
        assert a.dtype == original.dtype and a.shape == original.shape and np.array_equal(a, original)
        assert a.tobytes() == original.tobytes()
        assert not a.flags.writeable and not np.shares_memory(a, b)
        with pytest.raises(ValueError, match="read-only"):
            a[(0,) * a.ndim] = a[(0,) * a.ndim]
    assert type(first[3]) is np.float64 and first[3] == 1.25


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_validation.py::test_m2_7_arrays_keep_dtype_shape_and_bits_and_load_fresh_and_read_only')
def test_m2_8_a_dataclass_is_rebuilt_through_its_constructor_and_only_if_declared(package):
    """M2.8: a frozen dataclass payload is rebuilt by its constructor (its laws re-run: a payload whose field
    violates ``Spec``'s law is refused with ``Spec``'s own message), and only a type the function DECLARES it
    returns may be named by a payload (another class in the payload is refused, never imported)."""
    alpha = package.module("alpha")
    tree, arrays = api.encode_payload(alpha.Result(1.5, np.ones(3), 3, True))
    allowed = {f"{alpha.__name__}:Result"}
    rebuilt = api.decode_payload(tree, arrays, allowed)
    assert type(rebuilt) is alpha.Result and rebuilt.k == 1.5 and not rebuilt.field.flags.writeable
    with pytest.raises(ValueError, match="does not return"):
        api.decode_payload(tree, arrays, set())
    tree_spec, arrays_spec = api.encode_payload(alpha.Spec(2))
    tree_spec["fields"]["n"] = {"i": "-1"}
    with pytest.raises(ValueError, match="n is a count"):
        api.decode_payload(tree_spec, arrays_spec, {f"{alpha.__name__}:Spec"})


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_validation.py::test_m2_6_floats_cross_the_payload_bit_for_bit')
def test_m2_9_what_has_no_exact_payload_is_refused():
    """M2.9: an object array, a set, a list, a complex scalar, a function and a plain object have no exact,
    pickle-free payload: ``Unencodable``, never a lossy write."""
    unencodable = api.name("Unencodable")
    for value in (np.array([object()]), {1, 2}, [1.0], 1j, len, object()):
        with pytest.raises(unencodable):
            api.encode_payload(value)
