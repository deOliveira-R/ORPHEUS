r"""Step 4 of #405 P3: the clients — the two multi-region solves and the trajectory reading; the exact medium is not one.

Spec ``.claude/plans/reference_p3_spec.md`` §1.4, gates M4.1–M4.10. Client names are resolved through
``_traced_memo_api.CLIENTS``. The real-tree mutation witnesses (M4.4) run on a COPY of ``orpheus/``'s Python
files in ``tmp_path``, in a fresh interpreter whose ``sys.path`` starts at the copy; no tracked file is edited.
"""
from __future__ import annotations

import importlib
import json
import os
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from . import _traced_memo_api as api

pytestmark = pytest.mark.foundation

REPO = Path(__file__).resolve().parents[3]
ABA_DIR = REPO / "tests" / "gates" / "sn" / "verification" / "analytical"
SPHERE_Q = {"n_r": 8, "n_mu": 8, "n_traj_quad": 16}
SETTINGS = {"max_iter": 500, "tol": 1e-10, "initial_k": 1.0}


def _aba() -> Any:
    sys.path.insert(0, str(ABA_DIR))
    try:
        return importlib.import_module("_aba_reference")
    finally:
        sys.path.remove(str(ABA_DIR))


@pytest.fixture(scope="module")
def warm_root(tmp_path_factory):
    """One cache root for the module's cheap rows: the sphere and the exact medium are generated once."""
    return tmp_path_factory.mktemp("memo_cache")


def _sphere_reference():
    from orpheus.derivations.continuous.trajectory_resolvent.reference import trajectory_resolvent_reference
    from orpheus.geometry import CoordSystem

    return trajectory_resolvent_reference(_aba().aba_specification(CoordSystem.SPHERICAL), SPHERE_Q, **SETTINGS)


def _direct_sphere_kwargs(**overrides):
    aba = _aba()
    xs = aba.aba_xs_2g()
    return dict(radii=list(aba.ABA_RADII), sigma_t=xs.sigma_t, sigma_s=xs.sigma_s, nu_sigma_f=xs.nu_sigma_f, chi=xs.chi,
                **SPHERE_Q, **SETTINGS) | overrides


def _exact_specification():
    from orpheus.data.cells import CellCoefficient, Channel
    from orpheus.derivations.common.xs_library import get_mixture
    from orpheus.numerics.question import Eigen
    from orpheus.specification.specification import InfiniteMediumSpecification

    return InfiniteMediumSpecification(0, get_mixture("A", "2g"), Eigen(CellCoefficient.every(Channel.FISSION_EMISSION)))


def test_m4_1_the_clients_are_traced_memos():
    """M4.1 (a): the two multi-region solvers and the trajectory reading are bound as traced memos where every
    caller reaches them (the module attribute ``api.SOLVERS`` patches, and ``Billiard`` calls; the derivation's
    ``evaluate``, which ``ReferenceSolution.read`` calls)."""
    memo_type = api.name("TracedMemo")
    for client in api.CLIENTS:
        assert isinstance(api.client(client), memo_type), f"{client} is not a traced memo"


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_1_a_miss_generates_once_and_a_hit_never', 'tests/gates/numerics/test_traced_memo_clients.py::test_m4_1_the_clients_are_traced_memos')
def test_m4_1b_a_reference_reads_through_its_reading_memo(warm_root, monkeypatch):
    """M4.1 (b), the route: a sphere reference's first ``read`` writes one reading entry and one solve entry,
    starting interpreters; a SECOND reference over the same content reads with 0 interpreters started (the
    activation leg: the first count is not zero)."""
    from orpheus.numerics.observable import Eigenvalue

    spawns = api.SpawnCounter(monkeypatch)
    with api.cache_root(warm_root):
        _sphere_reference().read(Eigenvalue())
        cold = spawns.count
        _sphere_reference().read(Eigenvalue())
    assert cold >= 1, "the activation leg: the cold read generated nothing"
    assert spawns.count == cold, f"the warm read started {spawns.count - cold} interpreters"
    names = {api.entry_function(f) for f in api.entry_files(warm_root)}
    assert api.function_id("trajectory_reading") in names and api.function_id("solve_sphere") in names, names


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_3_the_signature_binds_one_call_to_one_key', 'tests/gates/numerics/test_traced_memo_process.py::test_m3_5_a_declared_canonical_form_is_both_the_key_and_what_the_child_receives')
def test_m4_2_one_solve_entry_serves_billiard_and_a_direct_caller(warm_root, monkeypatch):
    """M4.2: after the reference has read, a DIRECT call of the sphere solver, with the census's spellings (a
    list of radii, every array the caller's own, defaults omitted where ``Billiard`` passes them), is a ``Hit``
    and starts no interpreter; its k is the reference's k bit for bit (census Q1c: 1 key, 2 consumers)."""
    from orpheus.numerics.observable import Eigenvalue

    with api.cache_root(warm_root):
        k_reference = _sphere_reference().read(Eigenvalue()).value
        solver = api.client("solve_sphere")
        assert api.verdict_kind(solver.lookup(**_direct_sphere_kwargs())) == "Hit"
        spawns = api.SpawnCounter(monkeypatch)
        result = solver(**_direct_sphere_kwargs())
    assert spawns.count == 0
    assert float(result.k_eff).hex() == float(k_reference).hex()


def _reading_manifest(root: Path) -> dict:
    entries = [f for f in api.entry_files(root) if api.entry_function(f) == api.function_id("trajectory_reading")]
    assert entries, "no reading entry"
    return api.read_entry(entries[0])


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_9_arguments_cross_through_their_constructor', 'tests/gates/numerics/test_traced_memo_process.py::test_m3_12a_a_parent_entry_pins_its_child_by_reference')
def test_m4_3_the_reading_manifest_holds_the_construction_and_pins_the_solve_by_reference(warm_root):
    """M4.3, finding F1 on the real client: the reading entry's manifest holds the derivation's construction
    (``TrajectoryResolventDerivation.__post_init__``, ``Billiard.__post_init__``, ``_route``,
    ``_layered_xs_payload``, ``_read_isotropically``, ``reference_body``) and the rays' rebuild
    (``_SphereRays.per_group``), which a pickled derivation would have hidden (0 of 6 in the census's 2 of 2
    traced workloads); the solver's def is NOT in it (it ran in the solve child) and the solve child IS
    referenced."""
    from orpheus.numerics.observable import Eigenvalue

    with api.cache_root(warm_root):
        _sphere_reference().read(Eigenvalue())
    entry = _reading_manifest(warm_root)
    functions = {q for _, q, *_ in api.manifest_rows(entry, "DefPin")}
    construction = {"TrajectoryResolventDerivation.__post_init__", "Billiard.__post_init__", "_route",
                    "_layered_xs_payload", "_read_isotropically", "reference_body", "_SphereRays.per_group"}
    assert construction <= functions, sorted(construction - functions)
    assert "solve_greens_function_sphere_mr" not in functions
    assert [c[0] for c in api.manifest_rows(entry, "ChildPin")] == [api.function_id("solve_sphere")]


# ── M4.4: the real-tree witnesses, on a copy ─────────────────────────────────────

DRIVER = textwrap.dedent('''
    import json, sys
    copy, cache, action, aba_dir = sys.argv[1:5]
    sys.path[:0] = [copy, aba_dir]
    import orpheus
    assert orpheus.__file__.startswith(copy), orpheus.__file__
    import _aba_reference as aba
    from orpheus.geometry import CoordSystem
    from orpheus.numerics.observable import Eigenvalue
    from orpheus.derivations.continuous.trajectory_resolvent.reference import trajectory_resolvent_reference
    from orpheus.derivations.continuous.trajectory_resolvent.greens_function import solve_greens_function_sphere_mr
    from orpheus.numerics.traced_memo import cache_root
    Q = {"n_r": 8, "n_mu": 8, "n_traj_quad": 16}
    S = {"max_iter": 500, "tol": 1e-10, "initial_k": 1.0}
    xs = aba.aba_xs_2g()
    direct = dict(radii=list(aba.ABA_RADII), sigma_t=xs.sigma_t, sigma_s=xs.sigma_s, nu_sigma_f=xs.nu_sigma_f, chi=xs.chi, **Q, **S)
    with cache_root(cache):
        ref = trajectory_resolvent_reference(aba.aba_specification(CoordSystem.SPHERICAL), Q, **S)
        if action == "generate":
            print("K", ref.read(Eigenvalue()).value.hex())
        d = ref.derivation
        verdicts = {
            "reading": type(type(d).evaluate.lookup(d, Eigenvalue())).__name__,
            "solve": type(solve_greens_function_sphere_mr.lookup(**direct)).__name__,
        }
        print("VERDICTS", json.dumps(verdicts))
''')

#: (row, file under orpheus/, the def or statement to find, its neutral replacement, reading verdict, solve verdict)
TREE_WITNESSES = [
    ("construction-helper", "derivations/continuous/trajectory_resolvent/billiard.py",
     "def _layered_xs_payload(", None, "Stale", "Hit"),
    ("the-solver", "derivations/continuous/trajectory_resolvent/greens_function.py",
     "def solve_greens_function_sphere_mr(", None, "Stale", "Stale"),
    ("rays-of-the-other-geometry", "derivations/continuous/trajectory_resolvent/reference.py",
     "class _CylinderRays", "per_group", "Hit", "Hit"),
    ("reading-module-constant", "derivations/continuous/trajectory_resolvent/reference.py",
     '_RAYS: dict[str, type[_SphereRays] | type[_CylinderRays]] = {', "constant", "Stale", "Hit"),
]


def _neutral_edit(path: Path, anchor: str, method: str | None) -> None:
    """A behaviour-neutral AST edit: ``_ = None`` as the first statement of the def at ``anchor`` (or of the
    method ``method`` of the class at ``anchor``), or a trailing ``# `` -free duplicate key for a constant."""
    text = path.read_text()
    if text.count(anchor) != 1:
        raise AssertionError(f"{anchor!r} occurs {text.count(anchor)} times in {path}")
    if method == "constant":
        new = text.replace(anchor, anchor + '\n    "sphere_mr_alias": _SphereRays,')
    else:
        start = text.index(anchor)
        if method is not None:
            start = text.index(f"def {method}(", start)
        body = text.index(":\n", text.index(")", start)) + 2
        indent = len(text[body:]) - len(text[body:].lstrip(" "))
        if text[body:].lstrip(" ").startswith(('"""', "r\"\"\"")):  # step past the docstring
            close = text.index('"""', text.index('"""', body) + 3) + 4
            body = close
        new = text[:body] + " " * indent + "_ = None\n" + text[body:]
    path.write_text(new)


def _copy_tree(tmp_path: Path) -> Path:
    copy = tmp_path / "copy"
    for source in (REPO / "orpheus").rglob("*.py"):
        target = copy / source.relative_to(REPO)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    return copy


def _drive(copy: Path, cache: Path, action: str) -> dict:
    run = subprocess.run(
        [sys.executable, "-O", "-c", DRIVER, str(copy), str(cache), action, str(ABA_DIR)],
        capture_output=True, text=True, cwd=str(copy.parent), env={k: v for k, v in os.environ.items() if k != "PYTHONPATH"},
    )
    lines = [line for line in run.stdout.splitlines() if line.startswith("VERDICTS ")]
    if not lines:
        raise AssertionError(f"the driver failed:\n{run.stderr[-3000:]}")
    return json.loads(lines[-1][len("VERDICTS "):])


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_2_an_edit_to_what_ran_regenerates_and_nothing_else_does', 'tests/gates/numerics/test_traced_memo_process.py::test_m3_12b_a_mutation_only_a_child_traced_makes_the_parent_stale', 'tests/gates/numerics/test_traced_memo_clients.py::test_m4_3_the_reading_manifest_holds_the_construction_and_pins_the_solve_by_reference')
@pytest.mark.parametrize("row,file,anchor,method,reading,solve", TREE_WITNESSES, ids=[w[0] for w in TREE_WITNESSES])
def test_m4_4_an_edit_to_the_copied_tree_misses_exactly_where_the_reading_ran(tmp_path, row, file, anchor, method, reading, solve):
    """M4.4, the X1 witnesses on the REAL client, on a copy of ``orpheus/``: a neutral edit to a construction
    helper the reading process ran makes the reading stale and leaves the solve child valid; an edit to the
    solver makes the solve child stale and the reading stale THROUGH its reference (recursive validation); an
    edit to the cylinder's rays, which a sphere reading never ran, leaves both valid; a constant in the reading
    module's skeleton makes the reading stale."""
    copy = _copy_tree(tmp_path)
    cache = tmp_path / "cache"
    before = _drive(copy, cache, "generate")
    assert before == {"reading": "Hit", "solve": "Hit"}, before
    _neutral_edit(copy / "orpheus" / file, anchor, method)
    after = _drive(copy, cache, "lookup")
    assert after == {"reading": reading, "solve": solve}, f"{row}: {after}"


# ── the readings' values ────────────────────────────────────────────────────────


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_validation.py::test_m2_6_floats_cross_the_payload_bit_for_bit', 'tests/gates/numerics/test_traced_memo_clients.py::test_m4_1b_a_reference_reads_through_its_reading_memo')
def test_m4_5_a_served_reading_is_the_fresh_reading_bit_for_bit(warm_root):
    """M4.5: the memo's readings equal the in-process readings under ``bypass()`` bit for bit: the sphere's k,
    a point value and a flux integral of the fission production (``float.hex``)."""
    from orpheus.numerics.mesh_free_function import RegionwiseConstant
    from orpheus.numerics.observable import Eigenvalue, FluxIntegral, PointValue

    weight = FluxIntegral(RegionwiseConstant(_aba().aba_xs_2g().nu_sigma_f))
    observables = (Eigenvalue(), PointValue(position=1.0, group=1), weight)
    with api.cache_root(warm_root):
        served = [_sphere_reference().read(o).value.hex() for o in observables]
    with api.bypass():
        fresh = [_sphere_reference().read(o).value.hex() for o in observables]
    assert served == fresh


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_validation.py::test_m2_7_arrays_keep_dtype_shape_and_bits_and_load_fresh_and_read_only')
def test_m4_7_a_consumer_cannot_write_into_a_served_solve(warm_root):
    """M4.7, the positive control of the read-only acceptance run: the arrays of a served solve are read-only,
    so a consumer that writes into one raises where it would otherwise corrupt every later hit."""
    with api.cache_root(warm_root):
        result = api.client("solve_sphere")(**_direct_sphere_kwargs())
    for field in ("psi_g", "phi_g", "r_nodes", "mu_nodes", "region_at_node", "last_emission_density"):
        assert not getattr(result, field).flags.writeable, field
    with pytest.raises(ValueError, match="read-only"):
        result.phi_g[0, 0] = 0.0


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_11_a_raising_generator_writes_nothing_and_raises_its_own_error')
def test_m4_8_an_unconverged_solve_is_a_cached_value_and_its_reading_a_refusal_never_cached(tmp_path, monkeypatch):
    """M4.8: with ``max_iter=2`` the solve does not converge. The solve child is cached (``converged`` False is
    a value, read by 3 consumer sites); the reading's refusal crosses back as ``RuntimeError`` "did not
    converge" and writes nothing, so a second reading starts one interpreter (the reading) and not two."""
    from orpheus.derivations.continuous.trajectory_resolvent.reference import trajectory_resolvent_reference
    from orpheus.geometry import CoordSystem
    from orpheus.numerics.observable import Eigenvalue

    spec = _aba().aba_specification(CoordSystem.SPHERICAL)
    with api.cache_root(tmp_path):
        with pytest.raises(RuntimeError, match="did not converge"):
            trajectory_resolvent_reference(spec, SPHERE_Q, max_iter=2, tol=1e-10, initial_k=1.0).read(Eigenvalue())
        solve = api.client("solve_sphere")(**_direct_sphere_kwargs(max_iter=2))
        assert solve.converged is False and solve.iterations == 2
        spawns = api.SpawnCounter(monkeypatch)
        with pytest.raises(RuntimeError, match="did not converge"):
            trajectory_resolvent_reference(spec, SPHERE_Q, max_iter=2, tol=1e-10, initial_k=1.0).read(Eigenvalue())
    assert spawns.count == 1
    names = [api.entry_function(f) for f in api.entry_files(tmp_path)]
    assert names == [api.function_id("solve_sphere")], names


def test_m4_9_the_exact_medium_is_read_in_the_asking_process(tmp_path, monkeypatch):
    """M4.9, re-posed by the user's ruling of 2026-10-04 (Q5: "It can execute if it's fast"): the exact infinite
    medium reads in 0.01 s, against 0.11 s warm and 6.3 s cold through a generating process, so it is not
    memoised. Its construction and every reading start no interpreter and write no entry, its ``evaluate`` is a
    plain method, and its refusal of a point value is raised in this process. Its certificate is unaffected
    either way: ``ReferenceSolution`` re-checks every claim against the evaluation at construction."""
    from orpheus.derivations.common.exact_homogeneous import ExactInfiniteMediumDerivation, exact_infinite_medium_reference
    from orpheus.numerics.observable import Eigenvalue, PointValue

    assert not isinstance(ExactInfiniteMediumDerivation.__dict__["evaluate"], api.name("TracedMemo"))
    spawns = api.SpawnCounter(monkeypatch)
    with api.cache_root(tmp_path):
        reference = exact_infinite_medium_reference(_exact_specification())
        reference.read(Eigenvalue())
        with pytest.raises(ValueError, match="no position"):
            reference.derivation.evaluate(PointValue(position=0.5, group=0))
    assert spawns.count == 0
    assert not api.entry_files(tmp_path)


def test_m4_10_the_default_store_is_gitignored():
    """M4.10: the default cache root is ``.cache/references`` under the checkout, and git ignores it (a
    ``git add`` of the tree can never commit an entry). The control: a tracked path is not ignored."""
    root = api.name("default_root")()
    assert root == REPO / ".cache" / "references", root
    probe = subprocess.run(["git", "check-ignore", "-q", str(root / "x" / "key.npz")], cwd=REPO)
    control = subprocess.run(["git", "check-ignore", "-q", str(REPO / "orpheus" / "__init__.py")], cwd=REPO)
    assert probe.returncode == 0 and control.returncode == 1


# ── the characteristic reading: the successors of M4.1b-M4.8 ─────────────────────
#
# P1 step (e1b) of the characteristic-reference campaign: step (e) deletes the trajectory-resolvent family, and with
# it the three old clients these rows were posed on. ``CharacteristicDerivation.evaluate`` is the one memo client
# that survives, at ONE level: it solves inside the reading process, so its entry pins no child (``ChildPin`` empty).
# The successors, row by row (``scratch/characteristic_architecture/p1_step_e/ta_e1b/README_consumers.md``):
#
# * M4.1b (the route, the activation leg): ``tests/gates/derivations/test_characteristic_reference.py::
#   test_c3_an_equal_derivation_reads_from_the_memo_and_a_new_resolution_does_not`` (a cold read starts an
#   interpreter, an equal derivation's read starts none and reads the same bits);
# * M4.2 (one child entry, two consumers) and M4.4's ``the-solver`` witness: the child-memo contract, on a synthetic
#   client: ``test_traced_memo_process.py::test_m3_12e_a_child_generated_by_its_parent_serves_a_direct_caller`` and
#   ``::test_m3_12b_a_mutation_only_a_child_traced_makes_the_parent_stale``;
# * M4.7 (served arrays read-only): ``test_traced_memo_process.py::test_m3_14_a_served_payload_is_fresh_read_only_
#   and_bit_identical`` (a characteristic reading serves a float, no array);
# * M4.3, M4.4's other three witnesses, M4.5 and M4.8: the rows below.

_CHARACTERISTIC_DRIVER_SETUP = '''
from orpheus.data.cells import CellCoefficient, Channel
from orpheus.data.materials import Materials
from orpheus.derivations.common.xs_library import get_xs, make_mixture
from orpheus.derivations.continuous.characteristic import CharacteristicDerivation, Resolution, TransportResolution
from orpheus.geometry import BC, StructuredGeometry
from orpheus.numerics.question import Eigen
from orpheus.specification.specification import GeometrySpecification


def tiny_sphere(radius=1.0):
    """A one-region sphere of mixture A (two groups, P0 scattering) under vacuum, at the smallest resolution."""
    xs = get_xs("A", "2g")
    mixture = make_mixture(sig_t=xs["sig_t"], sig_c=xs["sig_c"], sig_f=xs["sig_f"], nu=xs["nu"], chi=xs["chi"],
                           sig_s=xs["sig_s"])
    geometry = StructuredGeometry.sphere((0.0, radius), (0,), outer=BC.vacuum)
    question = Eigen(CellCoefficient.every(Channel.FISSION_EMISSION))
    specification = GeometrySpecification(Materials({0: mixture}), geometry, question)
    return CharacteristicDerivation(specification, Resolution(1, 0, 0.4, TransportResolution(4, 4, 4), 4))
'''

_tiny_namespace: dict[str, Any] = {}
exec(_CHARACTERISTIC_DRIVER_SETUP, _tiny_namespace)


def _tiny_characteristic_sphere(radius: float = 1.0) -> Any:
    """The characteristic sphere the rows below read (the same text the copied-tree driver execs)."""
    return _tiny_namespace["tiny_sphere"](radius)


def _characteristic_entries(root: Path) -> list[Path]:
    return [f for f in api.entry_files(root) if api.entry_function(f) == api.function_id("characteristic_reading")]


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_9_arguments_cross_through_their_constructor',
                      'tests/gates/derivations/test_characteristic_reference.py::test_c3_an_equal_derivation_reads_from_the_memo_and_a_new_resolution_does_not')
def test_m4_3c_the_characteristic_reading_manifest_holds_its_construction_and_no_child(tmp_path):
    """M4.3, re-posed on the characteristic reading: the entry's manifest holds the derivation's construction
    (``CharacteristicDerivation.__post_init__``, ``Walls.of``, ``PanelBasis.of``, ``RegionCrossSections.of``,
    ``GalerkinSystem.__post_init__``) and the solve that ran in the reading process (``LinePeriod.inflow``,
    ``DensePencil.fundamental``), none of which a pickled derivation would show; the cylinder's polar rule, which a
    sphere never runs, is NOT in it; and it pins no child (the reference has one memo level). First red: the
    derivation crossing to the child pickled (its construction then runs in the parent and the DefPins lose it)."""
    from orpheus.numerics.observable import Eigenvalue

    with api.cache_root(tmp_path):
        _tiny_characteristic_sphere().evaluate(Eigenvalue())
    (entry_file,) = _characteristic_entries(tmp_path)
    entry = api.read_entry(entry_file)
    functions = {q for _, q, *_ in api.manifest_rows(entry, "DefPin")}
    construction = {"CharacteristicDerivation.__post_init__", "Walls.of", "PanelBasis.of", "RegionCrossSections.of",
                    "GalerkinSystem.__post_init__", "LinePeriod.inflow", "DensePencil.fundamental"}
    assert construction <= functions, sorted(construction - functions)
    assert "polar_rule" not in functions and "impact_per_polar" not in functions
    assert api.manifest_rows(entry, "ChildPin") == []


_CHARACTERISTIC_DRIVER = textwrap.dedent('''
    import sys
    copy, cache, action = sys.argv[1:4]
    sys.path[:0] = [copy]
    import orpheus
    assert orpheus.__file__.startswith(copy), orpheus.__file__
''') + _CHARACTERISTIC_DRIVER_SETUP + textwrap.dedent('''
    from orpheus.numerics.observable import Eigenvalue
    from orpheus.numerics.traced_memo import cache_root
    with cache_root(cache):
        d = tiny_sphere()
        if action == "generate":
            print("K", d.evaluate(Eigenvalue()).value.hex())
        print("VERDICT", type(type(d).evaluate.lookup(d, Eigenvalue())).__name__)
''')

#: (row, file under orpheus/, the def or statement to find, the edit, the reading's verdict after it)
CHARACTERISTIC_TREE_WITNESSES = [
    ("construction-helper", "derivations/continuous/characteristic/walls.py", "def _wall_of(", None, "Stale"),
    ("rule-of-another-chart", "derivations/continuous/characteristic/lines.py", "def polar_rule(", None, "Hit"),
    ("reading-module-constant", "derivations/continuous/characteristic/reference.py",
     'Lift = Callable[["sympy.Expr"], "sympy.Expr"]', "constant-before", "Stale"),
]


def _drive_characteristic(copy: Path, cache: Path, action: str) -> str:
    run = subprocess.run(
        [sys.executable, "-O", "-c", _CHARACTERISTIC_DRIVER, str(copy), str(cache), action],
        capture_output=True, text=True, cwd=str(copy.parent), env={k: v for k, v in os.environ.items() if k != "PYTHONPATH"},
    )
    lines = [line for line in run.stdout.splitlines() if line.startswith("VERDICT ")]
    if not lines:
        raise AssertionError(f"the driver failed:\n{run.stderr[-3000:]}")
    return lines[-1][len("VERDICT "):]


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_2_an_edit_to_what_ran_regenerates_and_nothing_else_does',
                      'tests/gates/numerics/test_traced_memo_clients.py::test_m4_3c_the_characteristic_reading_manifest_holds_its_construction_and_no_child')
@pytest.mark.parametrize("row,file,anchor,method,reading", CHARACTERISTIC_TREE_WITNESSES,
                         ids=[w[0] for w in CHARACTERISTIC_TREE_WITNESSES])
def test_m4_4c_an_edit_to_the_copied_tree_misses_exactly_where_the_characteristic_reading_ran(tmp_path, row, file, anchor,
                                                                                            method, reading):
    """M4.4, re-posed on the characteristic reading, on a copy of ``orpheus/``: a neutral edit to a construction
    helper the reading ran (``_wall_of``) makes the reading stale; an edit to the cylinder's polar rule, which a
    sphere reading never runs, leaves it a hit; a new constant in the reading module's skeleton makes it stale.
    The copy's own driver reports ``Hit`` before the edit (the activation leg)."""
    copy = _copy_tree(tmp_path)
    cache = tmp_path / "cache"
    assert _drive_characteristic(copy, cache, "generate") == "Hit"
    path = copy / "orpheus" / file
    if method == "constant-before":
        text = path.read_text()
        if text.count(anchor) != 1:
            raise AssertionError(f"{anchor!r} occurs {text.count(anchor)} times in {path}")
        path.write_text(text.replace(anchor, "_NEUTRAL_CONSTANT = None\n" + anchor))
    else:
        _neutral_edit(path, anchor, method)
    assert _drive_characteristic(copy, cache, "lookup") == reading, row


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_validation.py::test_m2_6_floats_cross_the_payload_bit_for_bit',
                      'tests/gates/derivations/test_characteristic_reference.py::test_c3_an_equal_derivation_reads_from_the_memo_and_a_new_resolution_does_not')
def test_m4_5c_a_served_characteristic_reading_is_the_fresh_reading_bit_for_bit(tmp_path):
    """M4.5, re-posed: the memo's characteristic readings equal the in-process readings under ``bypass()`` bit for
    bit: k, a point value of group 1 and a flux integral of the fission production (``float.hex``). The second
    pass through the memo reads each from an entry (three entries written)."""
    from orpheus.numerics.mesh_free_function import RegionwiseConstant
    from orpheus.numerics.observable import Eigenvalue, FluxIntegral, PointValue

    production = _tiny_characteristic_sphere().specification.materials[0]
    weight = FluxIntegral(RegionwiseConstant(np.asarray(production.SigP)[None, :]))
    observables = (Eigenvalue(), PointValue(position=0.4, group=1), weight)
    with api.cache_root(tmp_path):
        generated = [_tiny_characteristic_sphere().evaluate(o).value.hex() for o in observables]
        served = [_tiny_characteristic_sphere().evaluate(o).value.hex() for o in observables]
    assert len(_characteristic_entries(tmp_path)) == 3
    with api.bypass():
        fresh = [_tiny_characteristic_sphere().evaluate(o).value.hex() for o in observables]
    assert served == generated == fresh


@pytest.mark.rests_on('tests/gates/numerics/test_traced_memo_process.py::test_m3_11_a_raising_generator_writes_nothing_and_raises_its_own_error',
                      'tests/gates/derivations/test_characteristic_reference.py::test_c4_a_refusal_crosses_the_memo_with_its_type')
def test_m4_8c_a_characteristic_refusal_is_never_cached(tmp_path, monkeypatch):
    """M4.8, re-posed: the characteristic reference is a direct solve (no iteration to leave unconverged), so its
    refusal after the solve is the one M4.8's contract applies to: the flux of a mode nearest tau. Read through the
    memo it crosses back as ``NotImplementedError`` with its message and writes nothing, so a second reading starts
    an interpreter again; the mode's eigenvalue, read beside it, is cached. First red: a refusal cached as a value
    (the second read starts no interpreter)."""
    from dataclasses import replace

    from orpheus.numerics.observable import Eigenvalue, PointValue
    from orpheus.numerics.question import Nearest

    base = _tiny_characteristic_sphere()
    question = replace(base.specification.question, mode=Nearest(tau=0.05))
    derivation = type(base)(replace(base.specification, question=question), base.resolution)
    with api.cache_root(tmp_path):
        with pytest.raises(NotImplementedError, match="a higher mode has no flux scale"):
            derivation.evaluate(PointValue(position=0.4, group=0))
        spawns = api.SpawnCounter(monkeypatch)
        with pytest.raises(NotImplementedError, match="a higher mode has no flux scale"):
            type(base)(derivation.specification, base.resolution).evaluate(PointValue(position=0.4, group=0))
        refused_again = spawns.count
        derivation.evaluate(Eigenvalue())
    assert refused_again >= 1, f"the second refusal started {refused_again} interpreters: it was served from an entry"
    assert len(_characteristic_entries(tmp_path)) == 1
