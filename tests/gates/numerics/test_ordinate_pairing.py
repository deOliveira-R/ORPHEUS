r"""Gates of :func:`orpheus.numerics.quadrature.ordinate_pairing.ordinate_pairing`
(unit L4 of workflow WF2.2, layer DAG wave 2; plan
``.claude/plans/wave2_boundary_carve.md`` section 2.5 and 5.4).

``ordinate_pairing(quadrature, motion) -> Permutation`` is the equivariance
certification of a quadrature under a rigid motion, keyed on
``(quadrature, motion)`` with no law argument and no boundary vocabulary
(ruling 7 of ``layer_dag.md``). The mirror plane is stated ONCE, by the motion:
its normal is the ``-1`` eigenvector of ``motion.linear``, and the
inflow-to-outflow crossing is checked on that axis.

It folds onto :meth:`Quadrature.ordinate_permutation` (``None`` becomes
``QuadratureNotClosedUnderMotionError``) and asserts on the RETURNED
permutation the three claims that keep their teeth independently of the
matcher's built-in certificate: ERR-042 (the trace measure ``w |mu_a|``,
relative tolerance 1e-12), ERR-044 (involution) and ERR-045 (inflow to
outflow on the named axis).

Ladder (``rests_on`` edges are declared on the tests)::

    test_a_quadrature_not_closed_under_the_motion_is_refused     foundation
    test_ordinate_pairing_returns_the_certified_permutation      edge: equal to the matcher
        ^ real-table positive control
    test_the_three_errors_keep_their_teeth                       the evidence (ERR-042/044/045)

Gate 1 is NOT the evidence: it is equal to the matcher by shared upstream
(X4). The evidence is gate 2: three tables, each wrong in exactly ONE
invariant, each reddening exactly its own error class, taken from
``tests/gates/geometry/test_bc_universal_invariants.py`` and MOVED into this file (wave 2, slice S2b, rider (a): the
geometry gate retires with the laws at WF2.4, and a gate must not import another test module).

The surrogate quadrature ``_TableQuad`` of this file serves the surface
``n_ordinates``, ``weights``, ``measure``, ``axis_cosines``,
``mean_axis_cosine`` and ``ordinate_permutation``; ``ordinate_pairing`` may
read nothing else from a quadrature (the certification reads
``mean_axis_cosine``, the flux-question spelling, for the trace measure).

Mutations for the reviewers (each on a COPY of ``ordinate_pairing.py``, built
by transforming ``inspect.getsource``, never by hand; the verdict is per arm).
Red sets are MEASURED (`[M]` battery of the L4 review, 2026-10-23):

* M-measure: delete the measure comparison. 2 reds:
  ``test_the_three_errors_keep_their_teeth[weight_class_mispairing]`` and its
  realize-time twin.
* M-involution: delete the ``perm[perm] == arange`` comparison. 1 red:
  ``[mirror_after_same_class_swap]`` (the realize-time twin carries the
  involution row too, since the amendment).
* M-crossing: delete the sign-class comparison. 3 reds: ``[self_map_identity]``,
  its realize-time twin and ``test_the_axis_is_the_one_the_motion_names``.
* M-axis: derive the crossing axis from a hard-coded ``"x"``. 11 reds: the
  closure test, 8 gate-1 y/z rows, ``test_the_axis_is_the_one_the_motion_names``
  and the positive control.
* M-none: return the permutation without raising on ``None``. 1 red:
  ``test_a_quadrature_not_closed_under_the_motion_is_refused``.
* M-copy (positive control, many reds): return ``np.arange`` as the
  permutation. 10 reds; the gl8 y and z rows stay green (the permutation there
  is the identity).
* M-shape-guard: delete the diagonal guard of ``_mirror_axis``. Reddens the
  ``rotoreflection`` row of
  ``test_a_motion_that_is_not_an_axis_reflection_is_refused`` (before the
  amendment no row separated the guard from nothing).
"""

from __future__ import annotations

from typing import cast

import numpy as np
import pytest

from orpheus.geometry.faces import AXIS_NAMES
from orpheus.geometry.transformation import Permutation, RigidMotion
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.tangential import TANGENTIAL_EPS

#: This file's own path: the target of the ``rests_on`` edges below.
_THIS_FILE = "tests/gates/numerics/test_ordinate_pairing.py"


def _pairing():
    """The function under test, imported lazily: at the base tip the module
    does not exist, and a collection error would hide which gate is red."""
    from orpheus.numerics.quadrature.ordinate_pairing import ordinate_pairing

    return ordinate_pairing


def _errors():
    import importlib

    # importlib, not `from orpheus.numerics.quadrature import ordinate_pairing`:
    # the package re-exports the FUNCTION under the module's name.
    return importlib.import_module("orpheus.numerics.quadrature.ordinate_pairing")


def _mirror(axis: str) -> RigidMotion:
    """The reflection whose plane is normal to the named axis (``mirror_motion``
    is unit L1's and is not on this branch)."""
    return RigidMotion.reflection(normal=np.eye(3)[AXIS_NAMES.index(axis)])


#: The four quadratures of the B3.4c measurement, lazily built.
_QUADRATURES = {
    "gl8": lambda: Quadrature.gauss_legendre(n_ordinates=8),
    "product(2,4)": lambda: Quadrature.product(n_mu=2, n_phi=4),
    "ls6": lambda: Quadrature.level_symmetric(sn_order=6),
    "lebedev17": lambda: Quadrature.lebedev(17),
}


# ─────────────────────────────────────────────────────────────────────
# The poison fixtures (moved here from ``tests/gates/geometry/test_bc_universal_invariants.py``, WF2.2 rider (a))
# ─────────────────────────────────────────────────────────────────────
#
# Gate 2 re-poses the three tables of the geometry gate VERBATIM. Until wave 2's slice S2b they were imported from that
# file; WF2.4 deletes it with the laws it tests, so the numerics gate must not depend on a file that is going to be
# retired (a test importing another test module is also a collection hazard). The bodies are the geometry file's, with
# one change of spelling: the true mirror is this file's ``_mirror(axis)`` (a ``RigidMotion.reflection``) where the
# original read ``SelfPairedDeck.mirror(...).motion`` (the same motion).


class _TableQuad:
    """Quadrature stand-in: REAL nodes/weights, injectable specular
    pairing.

    The negative legs below need pairings that are wrong in exactly ONE
    invariant while every other datum stays production-real — the
    surgical way to prove the three pairing invariants are
    independent (the ERR-045 catalog lesson: "checking only one or two
    leaves a hole"). Served through ``ordinate_permutation`` — the
    certification's read since G6.3 step 7d — so every mutant must be a
    genuine bijection (:class:`Permutation` refuses anything else at
    construction), which the three mutant families all are: each
    violates its ONE invariant while remaining a permutation.
    """

    def __init__(self, base: Quadrature, table: np.ndarray) -> None:
        self.n_ordinates = base.n_ordinates
        self.weights = base.weights
        self._base = base
        self._table = np.asarray(table)

    @property
    def measure(self):
        # The surrogate mutates the PERMUTATION only; its nodes are the real
        # rule's, so the measure passes straight through. Needed because the
        # reference embeds via ``quad.measure.support.ambient_points(...)`` — a duck-typed
        # stand-in owes every member of the contract it stands in for, not
        # only the ones the last reader happened to call (``plan-authoring``
        # §6b: a surrogate is a call site no symbol grep can return).
        return self._base.measure

    def axis_cosines(self, axis_index: int) -> np.ndarray:
        return self._base.axis_cosines(axis_index)

    def mean_axis_cosine(self, axis_index: int) -> np.ndarray:
        return self._base.mean_axis_cosine(axis_index)

    def ordinate_permutation(self, motion) -> Permutation:  # noqa: ARG002
        return Permutation(self._table)


def _non_involutive_but_otherwise_valid(
    quadrature: Quadrature, *, axis: str,
) -> np.ndarray:
    r"""A table that breaks ONLY the involution (ERR-044).

    Why not the obvious ``np.roll(arange(N), 1)``: an n-cycle rotation breaks
    all three invariants at once (it pairs unequal weight classes, so ERR-042
    fires FIRST in the aggregate and the row proves nothing about ERR-044's
    independence). Worse, **no pure 3-cycle can ever be the ERR-044-only
    mutant**: a cycle :math:`a \to b \to c \to a` needs
    :math:`\mathrm{sign}(b) = -\mathrm{sign}(a)`,
    :math:`\mathrm{sign}(c) = -\mathrm{sign}(b) = \mathrm{sign}(a)` and then
    :math:`\mathrm{sign}(a) = -\mathrm{sign}(c) = -\mathrm{sign}(a)` — a
    contradiction, so every odd cycle breaks ERR-045 too.

    The construction that works is :math:`\pi \circ \sigma`, where :math:`\pi`
    is the TRUE mirror and :math:`\sigma` swaps two ordinates of the SAME sign
    and the SAME measure :math:`w|\mu_a|`:

    * **measure** — :math:`m[\pi[\sigma[i]]] = m[\sigma[i]] = m[i]`, since
      :math:`\sigma` stays inside a measure class and :math:`\pi` preserves it;
    * **sign** — :math:`\mathrm{sign}(\pi[\sigma[i]]) = -\mathrm{sign}(i)`,
      since :math:`\sigma` preserves sign and :math:`\pi` flips it;
    * **involution** — BROKEN: for the swapped pair :math:`\{a, b\}`,
      :math:`T[T[a]] = \pi[\sigma[\pi[b]]] = \pi[\pi[b]] = b \neq a`
      (:math:`\sigma` does not touch :math:`\pi[b]`, which has the opposite
      sign).
    """
    true_pi = quadrature.ordinate_permutation(_mirror(axis))
    if true_pi is None:  # pragma: no cover — fixtures are mirror-closed
        raise AssertionError(f"{quadrature!r} has no {axis}-mirror pairing")
    perm = true_pi.indices
    # ``mean_axis_cosine``: the w|mu| cosine measure is a FLUX question, so a
    # suppressed axis contributes zero rather than having no answer. Mirrors
    # the production spelling in ``geometry/boundary/_specular.py``.
    mu = quadrature.mean_axis_cosine(AXIS_NAMES.index(axis))
    measure = quadrature.weights * np.abs(mu)
    # Group same-sign ordinates by measure; the first class with ≥2 members
    # supplies the swap. Rounded because the class is a float coincidence of
    # a symmetric node set, not an exact-compare question.
    classes: dict[tuple[int, float], list[int]] = {}
    for n in range(int(quadrature.n_ordinates)):
        classes.setdefault(
            (int(np.sign(mu[n])), round(float(measure[n]), 12)), []
        ).append(n)
    for members in classes.values():
        if len(members) >= 2:
            a, b = members[0], members[1]
            sigma = np.arange(int(quadrature.n_ordinates))
            sigma[a], sigma[b] = b, a
            return perm[sigma]
    raise AssertionError(
        f"{quadrature!r} has no same-sign measure class with two members, so "
        f"the ERR-044-only mutant is not constructible here — pick a "
        f"quadrature with a degenerate measure class."
    )


#: ``id -> (id, quadrature_factory, table_factory)``. Each table is production-real except for the ONE property it breaks;
#: the quadrature is per-row (the ERR-044-only mutant needs a degenerate measure class, and ``gauss_legendre(8)`` has none).
_BROKEN = {
    "weight_class_mispairing": ("weight_class_mispairing", lambda: Quadrature.gauss_legendre(n_ordinates=8), lambda q: np.array([1, 0, 3, 2, 5, 4, 7, 6])),
    "self_map_identity": ("self_map_identity", lambda: Quadrature.gauss_legendre(n_ordinates=8), lambda q: np.arange(q.n_ordinates)),
    "mirror_after_same_class_swap": ("mirror_after_same_class_swap", lambda: Quadrature.level_symmetric(sn_order=6), lambda q: _non_involutive_but_otherwise_valid(q, axis="x")),
}


# ─────────────────────────────────────────────────────────────────────
# Gate 3: no closure, no pairing (the foundation: the refusal of the matcher's None)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_a_quadrature_not_closed_under_the_motion_is_refused() -> None:
    r"""An odd-``n_phi`` product rule has no x-mirror closure.

    Red input (today's tree): ``product(2, 3).ordinate_permutation(x-mirror)``
    is ``None`` (``[M]`` 2026-10-23, the matcher finds no bijection), and the
    positive control on the same rule about y, where the closure exists,
    returns a ``Permutation``. The refusal is a SUBCLASS of
    ``OrdinatePairingError``, and the error text names the motion's axis.
    """
    quad = Quadrature.product(n_mu=2, n_phi=3)
    assert quad.ordinate_permutation(_mirror("x")) is None, (
        "premise: the odd-n_phi rule is not closed under the x-mirror"
    )
    assert quad.ordinate_permutation(_mirror("y")) is not None, (
        "positive control: the same rule IS closed under the y-mirror"
    )
    module = _errors()
    with pytest.raises(module.QuadratureNotClosedUnderMotionError) as info:
        _pairing()(quad, _mirror("x"))
    assert isinstance(info.value, module.OrdinatePairingError)
    assert isinstance(info.value, ValueError)
    # the rule closed under the y-mirror is certified, not refused
    _pairing()(quad, _mirror("y"))


# ─────────────────────────────────────────────────────────────────────
# Gate 1: the returned permutation is the matcher's, certified
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS_FILE}::test_a_quadrature_not_closed_under_the_motion_is_refused")
@pytest.mark.parametrize("axis", AXIS_NAMES)
@pytest.mark.parametrize("name", list(_QUADRATURES))
def test_ordinate_pairing_returns_the_certified_permutation(name: str, axis: str) -> None:
    r"""The returned permutation equals ``quadrature.ordinate_permutation(motion)``,
    is an involution, and maps inflow ordinates to outflow ordinates on the
    axis the motion names.

    Equal by shared upstream (it folds onto the matcher), so this is NOT the
    evidence of the three checks (gate 2 is); it pins that the function RETURNS
    the certified array (X3: one computation, no second derivation) and takes
    no axis argument. Red input: the function does not exist at the base tip.
    """
    quad = _QUADRATURES[name]()
    motion = _mirror(axis)
    expected = quad.ordinate_permutation(motion)
    assert expected is not None, f"{name} must be closed under the {axis}-mirror"
    got = _pairing()(quad, motion)
    assert isinstance(got, Permutation)
    np.testing.assert_array_equal(got.indices, expected.indices)
    n = int(quad.n_ordinates)
    np.testing.assert_array_equal(got.indices[got.indices], np.arange(n))
    mu = quad.mean_axis_cosine(AXIS_NAMES.index(axis))
    active = np.abs(mu) > TANGENTIAL_EPS
    if active.any():
        assert (np.sign(mu[got.indices][active]) == -np.sign(mu[active])).all()


@pytest.mark.foundation
def test_ordinate_pairing_takes_no_axis_argument() -> None:
    r"""``ordinate_pairing(quadrature, motion)`` has exactly two parameters
    (C5: the axis is stated once, by the motion). Red input: a signature with a
    third ``axis`` parameter, or no function at all."""
    import inspect

    params = list(inspect.signature(_pairing()).parameters)
    assert params == ["quadrature", "motion"], params


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS_FILE}::test_a_quadrature_not_closed_under_the_motion_is_refused")
def test_a_motion_that_is_not_an_axis_reflection_is_refused() -> None:
    r"""A rotation, a translation, the identity, an oblique reflection and an
    inversion each raise ``OrdinatePairingError`` (NOT one of the three
    certification subclasses, and not the closure error: the motion names no
    mirror plane). Positive control: the axis reflection on the same rule
    certifies. Red input: the function does not exist at the base tip; after
    the build, a version that reads the axis from ``-1`` eigenvalue count
    alone accepts the inversion (three ``-1`` eigenvalues) and reddens the
    ``inversion`` row."""
    module = _errors()
    quad = Quadrature.level_symmetric(sn_order=6)
    _pairing()(quad, _mirror("x"))  # positive control
    motions = {
        "rotation_z_90": RigidMotion.rotation_about_axis(axis=np.array([0.0, 0.0, 1.0]), angle=np.pi / 2),
        "translation": RigidMotion.translation_by(np.array([1.0, 0.0, 0.0])),
        "identity": RigidMotion.identity(3),
        "oblique_reflection": RigidMotion.reflection(normal=np.array([1.0, 1.0, 0.0])),
        "inversion": RigidMotion.inversion(3),
        # diag(-1, 1, 1) @ R_x(pi/2): linear diagonal (-1, 0, 0). One -1 on the
        # diagonal, yet no mirror: with the diagonal guard deleted it is read
        # as an x-mirror and raises a certification subclass (`[M]`).
        "rotoreflection": RigidMotion.reflection(normal=np.array([1.0, 0.0, 0.0]))
        @ RigidMotion.rotation_about_axis(axis=np.array([1.0, 0.0, 0.0]), angle=np.pi / 2),
    }
    for label, motion in motions.items():
        with pytest.raises(module.OrdinatePairingError) as info:
            _pairing()(quad, motion)
        assert not isinstance(
            info.value,
            (
                module.QuadratureNotClosedUnderMotionError,
                module.PairingNotMeasurePreservingError,
                module.PairingNotInvolutiveError,
                module.PairingDoesNotCrossError,
            ),
        ), f"{label}: a non-reflection must raise the base error, not a certification subclass"


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS_FILE}::test_ordinate_pairing_returns_the_certified_permutation")
def test_the_axis_is_the_one_the_motion_names() -> None:
    r"""The crossing is checked on the motion's axis, not on a default. The
    identity table is the ERR-045 poison on the x axis of GL-8 (every
    ordinate self-maps, ``mu_x != 0``); it is accepted on the y axis of the
    same rule, where every ordinate is tangential and the self-map is the
    geometrically correct image. Red input: an implementation that hard-codes
    the checked axis fails one of the two legs."""
    module = _errors()
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    mutant = _TableQuad(quad, np.arange(quad.n_ordinates))
    with pytest.raises(module.PairingDoesNotCrossError):
        _pairing()(cast(Quadrature, mutant), _mirror("x"))
    # Same table, the motion names y: every ordinate is tangential there.
    np.testing.assert_array_equal(
        _pairing()(cast(Quadrature, mutant), _mirror("y")).indices, np.arange(quad.n_ordinates)
    )


# ─────────────────────────────────────────────────────────────────────
# Gate 2: the three errors keep their teeth (the evidence)
# ─────────────────────────────────────────────────────────────────────

#: broken-table id (a key of ``_BROKEN`` in the geometry gate) -> (error class
#: name in the new module, the ERR entry).
_EXPECTED = {
    "weight_class_mispairing": ("PairingNotMeasurePreservingError", "ERR-042"),
    "self_map_identity": ("PairingDoesNotCrossError", "ERR-045"),
    "mirror_after_same_class_swap": ("PairingNotInvolutiveError", "ERR-044"),
}


def _broken_rows():
    return _BROKEN


@pytest.mark.l1
@pytest.mark.catches("ERR-042", "ERR-044", "ERR-045")
@pytest.mark.rests_on(f"{_THIS_FILE}::test_ordinate_pairing_returns_the_certified_permutation")
@pytest.mark.parametrize("table_id", list(_EXPECTED))
def test_the_three_errors_keep_their_teeth(table_id: str) -> None:
    r"""Each broken table (wrong in exactly ONE invariant: the fixtures of
    this file's ``_BROKEN``, taken from ``test_bc_universal_invariants.py``) raises EXACTLY its own
    error class, and none of the other two (``type(error) is expected``).

    Red input: the function does not exist at the base tip. After the build,
    deleting one check (the mutations of the module docstring) reddens exactly
    the row of that check.
    """
    module = _errors()
    row = _broken_rows()[table_id]
    quad = row[1]()
    table = row[2](quad)
    mutant = _TableQuad(quad, table)
    expected = getattr(module, _EXPECTED[table_id][0])
    with pytest.raises(module.OrdinatePairingError) as info:
        _pairing()(cast(Quadrature, mutant), _mirror("x"))
    assert type(info.value) is expected, (
        f"{table_id}: expected exactly {expected.__name__}, got {type(info.value).__name__}"
    )


@pytest.mark.l1
@pytest.mark.catches("ERR-042", "ERR-044", "ERR-045")
@pytest.mark.rests_on(f"{_THIS_FILE}::test_ordinate_pairing_returns_the_certified_permutation")
@pytest.mark.parametrize("table_id", ["weight_class_mispairing", "self_map_identity", "mirror_after_same_class_swap"])
def test_a_poisoned_quadrature_pairing_reddens_at_use_time(
    table_id: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    r"""The realize-time fixture of the geometry gate (a real ``Quadrature``
    whose ``ordinate_permutation`` is poisoned) reddens ``ordinate_pairing``
    itself with the same classes: the certification is a production read of the
    matcher's result, not a check on a surrogate. Red input: the function does
    not exist at the base tip."""
    module = _errors()
    row = _broken_rows()[table_id]
    quad = row[1]()
    table = row[2](quad)
    monkeypatch.setattr(quad, "ordinate_permutation", lambda motion: Permutation(table))
    expected = getattr(module, _EXPECTED[table_id][0])
    with pytest.raises(expected) as info:
        _pairing()(quad, _mirror("x"))
    assert type(info.value) is expected


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS_FILE}::test_ordinate_pairing_returns_the_certified_permutation")
def test_the_independence_triad_has_a_positive_control_per_quadrature() -> None:
    r"""Anti-pattern 11: the three negatives above prove the RAISING, so the
    real tables of the same three quadratures certify on all three axes (a
    certification that raised on every table would pass every negative row).
    Red input: the function does not exist at the base tip."""
    for table_id in _EXPECTED:
        quad = _broken_rows()[table_id][1]()
        for axis in AXIS_NAMES:
            if quad.ordinate_permutation(_mirror(axis)) is not None:
                _pairing()(quad, _mirror(axis))


@pytest.mark.foundation
def test_the_package_exports_the_function_and_its_five_errors() -> None:
    r"""``orpheus.numerics.quadrature`` exports ``ordinate_pairing`` and the
    base error with its four subclasses, and the subclasses derive from the
    base. Red input: nothing is exported at the base tip."""
    import orpheus.numerics.quadrature as package

    names = [
        "ordinate_pairing",
        "OrdinatePairingError",
        "QuadratureNotClosedUnderMotionError",
        "PairingNotMeasurePreservingError",
        "PairingNotInvolutiveError",
        "PairingDoesNotCrossError",
    ]
    for name in names:
        assert name in package.__all__, name
    for name in names[2:]:
        assert issubclass(getattr(package, name), package.OrdinatePairingError)


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS_FILE}::test_the_three_errors_keep_their_teeth")
def test_a_folded_rule_is_diagnosed_as_having_spent_the_mirror() -> None:
    r"""On a rule whose orbit space SPENT the mirror the pairing is the identity
    and the outflow half is not stored: the refusal names the orbit space, not
    a wrong pairing (R2 of #434, re-posed from ``test_invariance.py``
    ``TestR2Err045NamesTheOrbitSpace`` ahead of that gate's retirement at S2).
    Control: the unfolded sibling certifies. Red input: a crossing check that
    does not read ``quadrature.measure`` (the predecessor-free message says
    "wrong pairing" and omits "SPENT the mirror")."""
    module = _errors()
    with pytest.raises(module.PairingDoesNotCrossError) as info:
        _pairing()(Quadrature.folded_product(4, 8), _mirror("y"))
    message = str(info.value)
    assert "SPENT the mirror" in message, message
    assert "wrong pairing" not in message, message
    _pairing()(Quadrature.product(4, 8), _mirror("y"))
