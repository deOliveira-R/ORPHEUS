"""L4 corroboration: the characteristic reference against the trajectory resolvent it replaces, during the migration.

**What this is.** Code-to-code agreement between two references
(``vv-principles``: L4, no correctness content). The characteristic
reference's correctness is its own ladder (``test_characteristic_*.py``,
against closed forms, mpmath routes and Garcia's table). This file's value is
migration safety (P1 step (c) of ``.claude/plans/characteristic_reference_architecture.md``,
spec ``scratch/characteristic_architecture/p1_verification_spec.md`` §7 rows
G1-G5): before the SN rows of step (d) move from the old reference to the new
one, the new reference must sit inside the old one's OWN error at each SN
row's problem, so that no SN tolerance recomputed on the new reference can
loosen. Each tolerance below is the old reference's own error estimate, never
widened to fit (the user's rulings of 2026-10-09).

**No level marker.** These are temporary rows: each is deleted, with the old
family it compares against, in step (e) (``retirement-audit`` D.14), and the
file is gone when its last row is. ``pyproject.toml`` defines l0-l3 and
foundation only, and an L4 row carries none of them (P0's
``tests/gates/geometry/test_kernel_corroboration.py`` is the pattern).

**Independence, asserted before every comparison** (``tests/gates/_corroboration.py``):

- statically, the old family's 20 modules and its three test helpers name nothing of the new package, by
  import, attribute chain or module-name string, and nothing their imports
  reach transitively is a new module;
- at run time, evaluating the old side executes no code object of the new
  package (counted through ``sys.monitoring``, the named entry points kept for
  the refusal's message). The old side runs under ``traced_memo.bypass()``: its evaluation is
  ``@traced_memo``, which on a miss runs in a fresh interpreter and on a hit
  runs nothing, so a spy in this process would read 0 either way
  (``test_the_runtime_leg_sees_a_memoised_call_only_under_bypass`` shows both).

When either leg fails, the row reds with "delete it in the migration commit".

**The new side** runs under ``bypass()`` too, so that an in-process mutation of
its reading reaches it (the first reds, ``scratch/characteristic_architecture/p1_step_c/ta_c/``).
Its resolution in each row is the one the explorer measured
(``p1_step_c/explorer_c.md`` §7): p = 5 for the sphere and slab rows, the door
default for the cylinder, the reading gates' ``_GARCIA_RES`` for Garcia.

**The old side** is read at its own fixture where the tolerance is that
fixture's error (G1-G3, G5), and at its FINEST rungs where the tolerance is a
step between finest rungs (G4: slab x4, sphere (36, 36, 96)), never at the
SN rows' coarser fixtures.
"""
from __future__ import annotations

import functools
import logging
import sys
import textwrap
from collections.abc import Iterator
from pathlib import Path

import numpy as np
import pytest

from orpheus.data.cells import CellCoefficient, Channel
from orpheus.data.materials import Materials
from orpheus.derivations.continuous.characteristic.grading import graded_ends
from orpheus.derivations.continuous.characteristic import (
    Resolution,
    TransportResolution,
    Walls,
    characteristic_reference,
)
from orpheus.geometry import BC, CoordSystem, StructuredGeometry
from orpheus.geometry.boundary import AlbedoBoundary, SpecularReturn
from orpheus.numerics.observable import Eigenvalue, Observable, PointValue
from orpheus.numerics.question import Eigen
from orpheus.numerics.traced_memo import bypass, cache_root, traced_memo
from orpheus.reference.solution import ReferenceSolution
from orpheus.specification.specification import GeometrySpecification
from tests.gates import _corroboration as corroboration
from tests.gates.derivations._trajectory_resolvent_ladders import (
    CYLINDER_3REG_REFERENCE_K,
    GARCIA_CASE1_RESOLVENT_STEP,
    SPHERE_3REG_REFERENCE_K,
    ceil_one_significant_figure,
    richardson_error,
    sphere_3reg_reference_ladder_estimate,
)
from tests.gates.derivations.test_characteristic_reading import _GARCIA, _GARCIA_RES
from tests.gates.derivations.test_peierls_greens_function_garcia2021 import (
    GARCIA_2021_CASE1_PHI,
    GARCIA_2021_CASE1_R,
    GARCIA_TO_VARIANT_ALPHA_FACTOR,
    _phi_at,
    solve_case1,
)
from tests.gates.sn.regression import _generate_snapshots as snapshots
from tests.gates.sn.verification.analytical import _aba_reference as aba
from tests.gates.sn.verification.analytical.test_partial_reflector_resolvent import _cross_sections

#: Each comparison logs its reading against its tolerance (``-o log_cli=true --log-cli-level=INFO`` prints them).
_LOG = logging.getLogger(__name__)

_OLD_PACKAGE = "orpheus.derivations.continuous.trajectory_resolvent"
#: The old family's size, ``[M]`` 2026-10-09 (the explorer's census, ``p1_step_c/explorer/imports.py``).
_OLD_MODULES = 20

#: The new side: the package, and the entry points every computation of it constructs or reads through (the door,
#: its evaluation, and the constructor of each object the door builds). A route reaching only a free function of
#: the package (``grading``) or a class not listed is not counted: the static legs are the guard there.
NEW = corroboration.NewSide(
    modules=("orpheus.derivations.continuous.characteristic",),
    entry_points=(
        "CharacteristicDerivation.__post_init__",
        "CharacteristicDerivation.evaluate",
        "GalerkinSystem.__post_init__",
        "GroupTransport.__post_init__",
        "LineRule.of",
        "PanelBasis.__post_init__",
        "RegionCrossSections.__post_init__",
        "Walls.of",
        "Walls.__post_init__",
    ),
    label="the new reference",
)

#: The new reference's resolutions (``p1_step_c/explorer_c.md`` §7): p = 5, the finest rung the explorer probed, where
#: k moved 4e-9 (slab) and 1.5e-8 (sphere, from p = 3) between rungs; and the door default for the cylinder, whose
#: p = 3 solve already costs about 21 min ``[M]`` (the plan's #587 entry).
_P5 = Resolution(5, 4, 0.4, TransportResolution(16, 20, 20), 12)
_DOOR_DEFAULT = Resolution(3, 2, 0.4, TransportResolution(8, 12, 12), 8)

_FISSION = Eigen(CellCoefficient.every(Channel.FISSION_EMISSION))
#: The old power iteration's tolerance on the A|B|A bodies (``_aba_reference.aba_reference_at``); the old value read
#: here is the tabulated rung when it is within ten times it.
_ABA_SOLVE_TOL = 1e-9


@pytest.fixture
def new_calls(monkeypatch: pytest.MonkeyPatch) -> Iterator[dict[str, int]]:
    """The runtime spy on the new reference: its entry points and every code object of the package."""
    with corroboration.spy(monkeypatch, NEW) as counts:
        yield counts


#: The old side: the family, and the test helpers the rows reach it through (qa 2, 2026-10-09): the A|B|A harness,
#: the Garcia solve and spline reading, and the ERR-094 rows' cross sections.
_OLD_HELPERS = (
    "tests.gates.sn.verification.analytical._aba_reference",
    "tests.gates.derivations.test_peierls_greens_function_garcia2021",
    "tests.gates.sn.verification.analytical.test_partial_reflector_resolvent",
)


def _assert_the_old_side_is_independent() -> None:
    """The static precondition on the old side: the per-module leg over the family and the helpers, and over every
    first-party module (``orpheus`` and ``tests``) they import."""
    old = (*corroboration.family(_OLD_PACKAGE), *_OLD_HELPERS)
    corroboration.assert_independent(old, NEW)
    corroboration.assert_closure_independent(old, NEW, ("orpheus", "tests"))


# ── the precondition and its positive controls ───────────────────────────────


def test_the_old_family_names_nothing_of_the_new_reference() -> None:
    """The static precondition over the old side: the family's 20 modules, the 3 test helpers the rows read it through,
    and every first-party module they import, each read by the per-module leg.

    ``[M]`` 2026-10-09: the family's closure is 153 first-party modules, 0 refused (qa ``closure_full_probe.py``). The
    family's size is asserted, so that an empty or partial walk cannot read as clean (``instrument-doctrine`` X1).
    """
    assert len(corroboration.family(_OLD_PACKAGE)) == _OLD_MODULES
    _assert_the_old_side_is_independent()


_SHIMS = {
    # name: (source, the leg that must refuse it: "direct", "closure" or None for a declared blindness)
    "_corroboration_shim_attribute": (
        "from orpheus.derivations import continuous\n"
        "def old():\n    return continuous.characteristic.Resolution\n", "direct"),
    "_corroboration_shim_string": (
        "import importlib\n"
        "def old():\n    return importlib.import_module('orpheus.derivations.continuous.characteristic')\n", "direct"),
    "_corroboration_shim_getattr": (
        "from orpheus.derivations import continuous\n"
        "def old():\n    return getattr(continuous, 'characteristic').Walls\n", "direct"),
    "_corroboration_shim_entry": (
        "def old():\n    return 'orpheus.derivations.continuous.characteristic:Walls'\n", "direct"),
    "_corroboration_shim_late": (
        "def old():\n    from orpheus.derivations.continuous.characteristic.walls import Walls\n    return Walls\n", "direct"),
    "_corroboration_shim_bridge": ("from orpheus.derivations.continuous.characteristic import Walls\n", "direct"),
    "_corroboration_shim_helper": ("from _corroboration_shim_bridge import Walls\n", "closure"),
    "_corroboration_shim_docstring": (
        '"""Replaced by :mod:`orpheus.derivations.continuous.characteristic`."""\n', None),
    "_corroboration_shim_computed": (
        "import importlib\n"
        "def old():\n    return importlib.import_module('orpheus.derivations.continuous.' + 'characteristic')\n", None),
}


def test_the_static_legs_see_each_shape_they_claim_and_not_a_docstring(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Positive controls for the static legs (X1), one per shape, and the negative legs beside them.

    Real modules: the new package's door imports its siblings RELATIVELY, and the reading gates import the package
    ABSOLUTELY; both are refused. Synthetic modules (written to ``tmp_path``): an attribute chain on an imported
    parent package and its literal ``getattr`` spelling, a module name in a string and in an entry-point string
    (``"module:Name"``), a late import inside a function, and a helper module between the old side and the new one
    (refused by the transitive leg alone, which applies the per-module leg to the helper; the per-module leg on the
    old module passes it, which is why the transitive leg exists). A docstring naming the package is not refused (prose imports nothing), and a module name
    COMPUTED at run time is not refused either: that blindness is declared in ``tests/gates/_corroboration.py`` and is
    the runtime leg's to cover.
    """
    for real in ("orpheus.derivations.continuous.characteristic.reference", "tests.gates.derivations.test_characteristic_reading"):
        with pytest.raises(AssertionError, match=r"now imports the new reference \(\{'import'.*this corroboration row compares the new reference"):
            corroboration.assert_independent([real], NEW)
    for name, (source, _) in _SHIMS.items():
        (tmp_path / f"{name}.py").write_text(textwrap.dedent(source))
    monkeypatch.syspath_prepend(str(tmp_path))
    first_party = ("orpheus", *_SHIMS)
    try:
        for name, (_, leg) in _SHIMS.items():
            if leg == "direct":
                shape = {"attribute": "attribute", "getattr": "attribute", "string": "string", "entry": "string"}.get(name.rsplit("_", 1)[1], "import")
                with pytest.raises(AssertionError, match=rf"\{{'{shape}'.*delete it in the migration commit"):
                    corroboration.assert_independent([name], NEW)
            elif leg == "closure":
                corroboration.assert_independent([name], NEW)
                with pytest.raises(AssertionError, match=r"the old family's imports reach the new reference: "
                                                         r"_corroboration_shim_bridge <- _corroboration_shim_helper now imports"):
                    corroboration.assert_closure_independent([name], NEW, first_party)
            else:
                corroboration.assert_independent([name], NEW)
                corroboration.assert_closure_independent([name], NEW, first_party)
    finally:
        for name in _SHIMS:
            sys.modules.pop(name, None)


@traced_memo
def _memoised_spelling_that_builds_the_new_walls(radius: float) -> float:
    """A stand-in "old spelling" behind a traced memo that reaches the new reference (``Walls.of``) with no import of it
    in the row: the shape the runtime leg exists for."""
    return float(len(Walls.of(StructuredGeometry.sphere((0.0, radius), (0,), outer=BC.vacuum)).walls))


def test_the_runtime_leg_sees_a_memoised_call_only_under_bypass(new_calls: dict[str, int], tmp_path: Path) -> None:
    """Positive control for the runtime leg (X1), and the reason the old side runs under ``bypass()``.

    A free function of the package (``grading.graded_ends``), which no entry point fronts, is refused by name: the
    leg counts every code object of the package (qa A2, green before the leg did). Under ``bypass()`` (what :func:`~tests.gates._corroboration.without` does) the memoised stand-in runs here and
    the spy counts its ``Walls.of``: refused. Called through the memo instead (a fresh store in ``tmp_path``, so a
    miss), the same call runs in a fresh interpreter and the spy reads 0 while the answer arrives: the vacuous
    reading the bypass removes. A call to the new reference outside the leg is counted and not charged, and a clean
    callable passes.
    """
    with pytest.raises(AssertionError, match=r"the old spelling called the new reference \(\{'grading\.graded_ends': 1\}\)"):
        corroboration.without(new_calls, NEW, graded_ends, 0.0, 1.0, True, True, 3, 0.4)
    with pytest.raises(AssertionError, match="the old spelling called the new reference .*'Walls.of': 1.*'walls.Walls.of': 1"):
        corroboration.without(new_calls, NEW, _memoised_spelling_that_builds_the_new_walls, 2.0)
    before = dict(new_calls)
    with cache_root(tmp_path):
        assert _memoised_spelling_that_builds_the_new_walls(2.0) == 1.0
    assert new_calls == before, "a fresh interpreter's calls reached this process's spy"
    Walls.of(StructuredGeometry.sphere((0.0, 1.0), (0,), outer=BC.vacuum))
    assert new_calls["Walls.of"] == before["Walls.of"] + 1
    assert corroboration.without(new_calls, NEW, lambda x: x + 1.0, 1.0) == 2.0


# ── the comparisons ──────────────────────────────────────────────────────────


@functools.cache
def _new(specification: GeometrySpecification, resolution: Resolution) -> ReferenceSolution:
    """The new reference, one object per (specification, resolution) in a session, so its solve is shared by the rows."""
    return characteristic_reference(specification, resolution)


def _read_new(specification: GeometrySpecification, resolution: Resolution, observable: Observable) -> float:
    with bypass():
        return float(_new(specification, resolution).read(observable).value)


def _old_aba(coord: CoordSystem, observable: Observable) -> float:
    """The old reference's reading on the A|B|A body at its fixture, from a FRESH reference (``aba_reference_at``).

    Not the session-cached ``aba_reference``: the old side carries no in-process state from one row to the next, so
    each row's runtime leg sees every call its old reading makes. A value cached by an earlier row whose window was
    refused would pass this row's leg vacuously (``[M]`` 2026-10-09, battery arm R: with the Garcia solve cached, the
    surface row stayed green behind the interior row's refusal).
    """
    return float(aba.aba_reference_at(coord, aba.ABA_REFERENCE_QUADRATURE[coord]).read(observable).value)


@pytest.mark.slow
def test_the_aba_sphere_k_agrees_with_the_old_reference_within_its_ladder_estimate(new_calls: dict[str, int]) -> None:
    """G1: the A|B|A sphere's k, new (p = 5) against old (36, 96, 64), within the old fixture's ladder estimate.

    Tolerance: ``sphere_3reg_reference_ladder_estimate()["k"]``, 3.22e-4 relative (the Richardson error of the radial
    step 36 -> 72 at second order plus the alternating angular step 96 -> 192, over the old table). The old value is
    first held to its tabulated rung (``SPHERE_3REG_REFERENCE_K[(36, 96)]``, within ten times its solve tolerance),
    so the reading is the one the estimate describes. ``[M]`` 2026-10-09 (explorer, ``g12_sphere.log``): 8.33e-5.
    Cost: about 190 s old, 10 s new.
    """
    _assert_the_old_side_is_independent()
    k_old = corroboration.without(new_calls, NEW, _old_aba, CoordSystem.SPHERICAL, Eigenvalue())
    np.testing.assert_allclose(k_old, SPHERE_3REG_REFERENCE_K[(36, 96)], rtol=10 * _ABA_SOLVE_TOL, atol=0.0)
    k_new = _read_new(aba.aba_specification(CoordSystem.SPHERICAL), _P5, Eigenvalue())
    tolerance = sphere_3reg_reference_ladder_estimate()["k"]
    _LOG.info("G1: k_old %r, k_new %r, relative %.3e, tolerance %.3e", k_old, k_new, abs(k_new - k_old) / k_old, tolerance)
    assert abs(k_new - k_old) / k_old <= tolerance, (k_new, k_old, tolerance)


@pytest.mark.slow
def test_the_aba_sphere_shape_agrees_with_the_old_reference_within_its_ladder_estimate(new_calls: dict[str, int]) -> None:
    """G2: the A|B|A sphere's shape, new (p = 5) against old (36, 96, 64), within the old fixture's ladder estimate.

    The shape is the 80 fission-gauged cell averages on the SN rows' 40-cell mesh
    (``_aba_reference.shape_observables``); the metric is the largest ``|new - old|`` over cells and groups, divided
    by the largest old average M. Tolerance: ``sphere_3reg_reference_ladder_estimate()["shape"]``, 1.37e-3, whose
    steps were measured on the retired nodal reading (``_trajectory_resolvent_ladders.py``, the
    ``SPHERE_3REG_REFERENCE_SHAPE_STEP`` comment): a caveat the plan carries for step (d). ``[M]`` 2026-10-09
    (explorer): 4.05e-5. Cost: about 290 s old (its own solve, then the 80 readings), 11 s new.
    """
    _assert_the_old_side_is_independent()
    observables = aba.shape_observables(snapshots._sphere_3region("2g", 40)["mesh"])
    assert len(observables) == 80

    def old_shape() -> np.ndarray:
        reference = aba.aba_reference_at(CoordSystem.SPHERICAL, aba.ABA_REFERENCE_QUADRATURE[CoordSystem.SPHERICAL])
        return np.array([float(reference.read(ratio).value) for _, _, ratio in observables])

    old = corroboration.without(new_calls, NEW, old_shape)
    new = np.array([_read_new(aba.aba_specification(CoordSystem.SPHERICAL), _P5, ratio) for _, _, ratio in observables])
    metric = float(np.max(np.abs(new - old)) / np.max(np.abs(old)))
    tolerance = sphere_3reg_reference_ladder_estimate()["shape"]
    _LOG.info("G2: M %r, max|new - old|/M %.3e, tolerance %.3e", float(np.max(np.abs(old))), metric, tolerance)
    assert metric <= tolerance, (metric, tolerance)


#: G3's tolerance: the old cylinder's step 32 -> 64 azimuthal nodes at 8 axial nodes, relative
#: (``_trajectory_resolvent_ladders.CYLINDER_3REG_REFERENCE_K`` comment: 7.8e-3, 1.9e-3, 1.4e-3 over 16 -> 128; not
#: monotone, so no bound). A comment's number, not a function's: the old family computes none.
_G3_AZIMUTHAL_STEP = 1.9e-3


@pytest.mark.slow
def test_the_aba_cylinder_k_lies_within_the_old_references_azimuthal_step(new_calls: dict[str, int]) -> None:
    """G3: the A|B|A cylinder's k, new (the door default) against old (24, 16, 32, 64), within 1.9e-3 relative.

    The old family has no bound on this body, only a non-monotone azimuthal ladder, so this row can REFUTE the new
    reference (a k outside the old one's largest measured step) and never certify it. The old value is first held to
    ``CYLINDER_3REG_REFERENCE_K[(24, 16, 32)]`` within ten times its solve tolerance. ``[M]`` (the plan, #587):
    5.63e-4. Cost, cold: about 13 min old, 21 min new; run once, detached (``p1_step_c/ta_c/``).
    """
    _assert_the_old_side_is_independent()
    k_old = corroboration.without(new_calls, NEW, _old_aba, CoordSystem.CYLINDRICAL, Eigenvalue())
    np.testing.assert_allclose(k_old, CYLINDER_3REG_REFERENCE_K[(24, 16, 32)], rtol=10 * _ABA_SOLVE_TOL, atol=0.0)
    k_new = _read_new(aba.aba_specification(CoordSystem.CYLINDRICAL), _DOOR_DEFAULT, Eigenvalue())
    _LOG.info("G3: k_old %r, k_new %r, relative %.3e, tolerance %.3e", k_old, k_new, abs(k_new - k_old) / k_old, _G3_AZIMUTHAL_STEP)
    assert abs(k_new - k_old) / k_old <= _G3_AZIMUTHAL_STEP, (k_new, k_old)


def _partial_reflector(geometry: StructuredGeometry) -> GeometrySpecification:
    """The ERR-094 rows' body (mixture A, two groups) posed for the new door: ``isotropic_mixture("A")``, since the
    library's mixture A carries a P1 moment the new reference refuses and the old solver never reads."""
    return GeometrySpecification(Materials({0: aba.isotropic_mixture("A")}), geometry, _FISSION)


def _assert_the_old_inputs_are_the_posed_mixture() -> tuple[np.ndarray, ...]:
    """The old solver's arrays (the SN rows' ``_cross_sections``) are, bit for bit, the posed mixture's: one problem."""
    sigma_t, sigma_s, nu_sigma_f, chi = _cross_sections()
    mixture = aba.isotropic_mixture("A")
    assert np.array_equal(np.asarray(mixture.SigT, dtype=float), sigma_t)
    assert len(mixture.SigS) == 1 and np.array_equal(mixture.SigS[0].toarray(), sigma_s)
    assert np.array_equal(np.asarray(mixture.SigP, dtype=float), nu_sigma_f)
    assert np.array_equal(np.asarray(mixture.chi, dtype=float), chi)
    assert all(m.nnz == 0 for m in mixture.Sig2), "an (n,2n) channel the old solver does not read"
    return sigma_t, sigma_s, nu_sigma_f, chi


#: G4 slab: the ceiling the user ruled (2026-10-09) for the tolerance this row computes from the old ladder; a
#: computed tolerance above it means the old ladder moved, and the row refuses rather than widen.
_G4_SLAB_RULED = 4e-5


@pytest.mark.slow
def test_the_partial_reflector_slab_k_agrees_with_the_old_finest_rung_within_its_richardson_error(new_calls: dict[str, int]) -> None:
    """G4, slab: L = 4 cm, albedos 0.3 | 0.7, mixture A, new (p = 5) against the old finest rung x4 = (64, 96, 128).

    Tolerance: the Richardson error of the old FINEST rung, computed here from the old ladder's own rungs x2, x3, x4
    (``x_m = (16 m, 24 m, 32 m)`` in ``(n_x, n_mu, n_traj_quad)``) at the order OBSERVED on them (qa 3):
    ``richardson_error(x4 - x3, 4/3, p)`` is x3's error, and x4's is that less the step; relative, rounded up by
    ``ceil_one_significant_figure``, and never above the ruled 4e-5. The observed order must lie within 0.2 of the
    scheme's 2 (the error model's premise: rungs in the asymptotic range). The last step alone (2.5e-5, the spec's first figure) understates x4's error 1.3x.
    ``[M]`` 2026-10-09 (explorer): the error is 3.2e-5, so 4e-5; the new value is 3.28e-5 from x4, and matches the
    ladder's extrapolated limit x4 + 2.7e-5 to about 4e-7. The arrays are asserted equal to the old solver's inputs.
    Cost: about 130 s old (x4 72 s), 11 s new.
    """
    _assert_the_old_side_is_independent()
    sigma_t, sigma_s, nu_sigma_f, chi = _assert_the_old_inputs_are_the_posed_mixture()
    from orpheus.derivations.continuous.trajectory_resolvent.greens_function_slab_asymmetric import (
        solve_greens_function_slab_asymmetric_mg,
    )

    def old_rung(m: int) -> float:
        result = solve_greens_function_slab_asymmetric_mg(
            4.0, sigma_t, sigma_s, nu_sigma_f, chi, alpha_left=0.3, alpha_right=0.7,
            n_x=16 * m, n_mu=24 * m, n_traj_quad=32 * m, max_iter=500, tol=1e-12,
        )
        if not result.converged:
            pytest.fail(f"the old slab rung x{m} did not converge")
        return float(result.k_eff)

    x2, x3, x4 = (corroboration.without(new_calls, NEW, old_rung, m) for m in (2, 3, 4))
    observed_ratio = (x4 - x3) / (x3 - x2)
    order = _order_from_step_ratio(observed_ratio, (2, 3, 4))
    assert abs(order - 2.0) <= 0.2, f"the old slab ladder is not second order (observed {order:.3f}): the error model fails"
    estimate = (richardson_error(x4 - x3, 4 / 3, order) - abs(x4 - x3)) / x4     # x4's error at the observed order
    tolerance = ceil_one_significant_figure(estimate)
    assert tolerance <= _G4_SLAB_RULED, f"the old ladder's error grew to {tolerance:.1e}: re-rule, never widen"
    geometry = StructuredGeometry.slab((0.0, 4.0), (0,), left=AlbedoBoundary(0.3, SpecularReturn("x")),
                                       right=AlbedoBoundary(0.7, SpecularReturn("x")))
    k_new = _read_new(_partial_reflector(geometry), _P5, Eigenvalue())
    _LOG.info("G4 slab: x2 %r, x3 %r, x4 %r, order %.4f, estimate %.4e, k_new %r, relative %.4e, tolerance %.1e",
              x2, x3, x4, order, estimate, k_new, abs(k_new - x4) / x4, tolerance)
    assert abs(k_new - x4) / x4 <= tolerance, (k_new, x4, tolerance)


def _order_from_step_ratio(ratio: float, m: tuple[int, int, int]) -> float:
    r"""The order p with :math:`(m_1^{-p} - m_2^{-p}) / (m_0^{-p} - m_1^{-p})` equal to ``ratio`` (bisection on
    [0.5, 8], where it decreases in p): rungs with error :math:`C m^{-p}`."""
    def model(p: float) -> float:
        return (m[1] ** -p - m[2] ** -p) / (m[0] ** -p - m[1] ** -p)

    lo, hi = 0.5, 8.0
    if not model(hi) <= ratio <= model(lo):
        return float("nan")
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if model(mid) > ratio else (lo, mid)
    return 0.5 * (lo + hi)


#: G4 sphere's tolerance: the old step (24, 24, 64) -> (36, 36, 96), 1.03e-4 relative, rounded to 1.0e-4 by the
#: spec (``test_partial_reflector_resolvent.py`` module docstring table).
_G4_SPHERE_STEP = 1.0e-4


@pytest.mark.slow
def test_the_partial_reflector_sphere_k_agrees_with_the_old_finest_rung_within_its_last_step(new_calls: dict[str, int]) -> None:
    """G4, sphere: R = 4 cm, albedo 0.7, mixture A, new (p = 5) against the old finest rung (36, 36, 96), within 1.0e-4.

    Tolerance: the old ladder's last step, (24, 24, 64) -> (36, 36, 96), 1.03e-4 relative (the SN module's docstring
    table). ``[M]`` 2026-10-09 (explorer): 4.13e-5. The arrays are asserted equal to the old solver's inputs. Cost:
    about 10 s old, 1 s new.
    """
    _assert_the_old_side_is_independent()
    sigma_t, sigma_s, nu_sigma_f, chi = _assert_the_old_inputs_are_the_posed_mixture()
    from orpheus.derivations.continuous.trajectory_resolvent.greens_function import solve_greens_function_sphere_mg

    def old() -> float:
        result = solve_greens_function_sphere_mg(4.0, sigma_t, sigma_s, nu_sigma_f, chi, alpha=0.7,
                                                 n_r=36, n_mu=36, n_traj_quad=96, max_iter=500, tol=1e-11)
        if not result.converged:
            pytest.fail("the old sphere rung (36, 36, 96) did not converge")
        return float(result.k_eff)

    k_old = corroboration.without(new_calls, NEW, old)
    k_new = _read_new(_partial_reflector(StructuredGeometry.sphere((0.0, 4.0), (0,), outer=AlbedoBoundary(0.7, SpecularReturn("x")))),
                      _P5, Eigenvalue())
    _LOG.info("G4 sphere: k_old %r, k_new %r, relative %.3e, tolerance %.1e", k_old, k_new, abs(k_new - k_old) / k_old, _G4_SPHERE_STEP)
    assert abs(k_new - k_old) / k_old <= _G4_SPHERE_STEP, (k_new, k_old)


def _old_garcia() -> tuple[float, ...]:
    """The old reference's Garcia Case 1 flux at Table 5's 15 radii, (n_r, n_mu) = (48, 24), the nodal spline reading.

    Solved in each row's own window, never cached across rows (``_old_aba``'s reason).
    """
    solution = solve_case1(48, 24)
    if not solution.converged:
        pytest.fail("the old Garcia solve did not converge")
    return tuple(_phi_at(solution, float(r)) for r in GARCIA_2021_CASE1_R)


def _new_garcia(radii: np.ndarray) -> np.ndarray:
    return np.array([_read_new(_GARCIA, _GARCIA_RES, PointValue(float(r), 0)) for r in radii])


def test_garcias_interior_fluxes_agree_with_the_old_reference_within_its_ladder_step(new_calls: dict[str, int]) -> None:
    """G5, interior: Garcia 2021 Case 1 at the 14 table radii inside the sphere, new (``_GARCIA_RES``) against old.

    Tolerance: ``GARCIA_CASE1_RESOLVENT_STEP["interior"]``, 1.43e-3 relative, the old fixture's largest change to its
    finer rungs over these radii. ``[M]`` 2026-10-09 (explorer): at most 1.07e-3, at r = 5.0 cm (an interface). Cost:
    about 6 s old, 6 s new.
    """
    _assert_the_old_side_is_independent()
    old = np.array(corroboration.without(new_calls, NEW, _old_garcia))
    interior = GARCIA_2021_CASE1_R < GARCIA_2021_CASE1_R[-1]
    assert int(interior.sum()) == 14
    new = _new_garcia(GARCIA_2021_CASE1_R[interior])
    relative = np.abs(new - old[interior]) / np.abs(old[interior])
    tolerance = GARCIA_CASE1_RESOLVENT_STEP["interior"]
    _LOG.info("G5 interior: max relative %.3e at r = %.1f, tolerance %.2e", relative.max(), GARCIA_2021_CASE1_R[interior][relative.argmax()], tolerance)
    assert relative.max() <= tolerance, dict(zip(GARCIA_2021_CASE1_R[interior].tolist(), relative.tolist()))


#: G5 surface: the ceiling the user ruled (2026-10-09) for the tolerance this row computes.
_G5_SURFACE_RULED = 3e-2


def test_garcias_surface_flux_agrees_with_the_old_reference_within_its_distance_from_garcia(new_calls: dict[str, int]) -> None:
    """G5, surface: Garcia 2021 Case 1 at r = R = 7 cm, new (``_GARCIA_RES``) against old.

    Tolerance: the old value's measured distance from Garcia's printed value there (half of it, the old reference's
    convention), relative, rounded up by ``ceil_one_significant_figure``: the old reference's own error at the
    surface, which its ladder step (1.65e-2, the spec's first figure) understates. ``[M]`` 2026-10-09 (explorer): the
    old value 2.085325 is 2.3e-2 from Garcia, so 3e-2; the new value 2.038175 is 1.2e-5 from Garcia and 2.26e-2 from
    the old one. The step (d) rows that rest on the old step are re-pointed there (the ruling).
    """
    _assert_the_old_side_is_independent()
    old = corroboration.without(new_calls, NEW, _old_garcia)[-1]
    garcia = GARCIA_TO_VARIANT_ALPHA_FACTOR * float(GARCIA_2021_CASE1_PHI[-1])
    tolerance = ceil_one_significant_figure(abs(old - garcia) / garcia)
    assert tolerance <= _G5_SURFACE_RULED, f"the old surface value moved {abs(old - garcia) / garcia:.1e} from Garcia: re-rule, never widen"
    new = float(_new_garcia(GARCIA_2021_CASE1_R[-1:])[0])
    _LOG.info("G5 surface: old %r, new %r, Garcia/2 %r, relative %.3e, tolerance %.0e", old, new, garcia, abs(new - old) / abs(old), tolerance)
    assert abs(new - old) / abs(old) <= tolerance, (new, old, tolerance)
