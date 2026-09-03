r"""#434 R2 — invariance is the measure's question: the gates.

Designed BEFORE the implementation (test-architect, 2026-09-03; the plan is
``scratch/_r2_verification_plan.md``, untracked) and landed with it.  EVERY
class lives here: the import graph (A), ONE closure (B), the deleted step 2
(C), the position and azimuth windows (D, E), candidates from the barycentres
(F), strict maximality (G), the retired names and the measure verbs (H), the
frozen ordinate-permutation and selection tables (H4), and ERR-045's fold
diagnosis (I).  The four entry points the old import cycle killed are a
RE-KEY of ``tests/test_layer_imports.py``'s fresh-interpreter list, not rows
here.  The mutation battery that proves these rows bite is
``scratch/_r2_mut.py`` (23 arms; its arm table lives there, once).

**Markers.**  Every row here is a SOFTWARE INVARIANT with no theory
``:label:``, so it is ``foundation`` and carries NO ``verifies``

**Discipline honoured here.**  No production call in any ``parametrize`` argument list
(rows parametrize by a LABEL and build in the body — `vv` Mode-8 third pipeline class);
every ``pytest.raises`` names a production entry point and a ``match=`` fragment; every
spy asserts its own installation; ``maximal_invariance_groups``/``symmetry_groups``
answers are compared as SETS, never as tuples (`[M]` walk and bruteforce return different
tuple ORDER on a chart-spelled fold).
"""

from __future__ import annotations

import ast
import pathlib
import subprocess
import sys

import numpy as np
import pytest

from orpheus.geometry.boundary import SelfPairedDeck
from orpheus.geometry.boundary._errors import (
    ReflectionDidNotMapInflowToOutflowError,
)
from orpheus.geometry.boundary._specular import (
    assert_specular_pairing_maps_inflow_to_outflow,
)
from orpheus.numerics import invariance as inv
from orpheus.numerics import manifold as man
from orpheus.numerics import symmetry as sym
from orpheus.numerics.invariance import (
    _distinct_azimuths,
    _embedded_nodes,
    _maximal,
    candidate_groups,
    certificate_under,
    symmetry_groups,
)
from orpheus.numerics.manifold import SPHERE, Quotient
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.quadrature.directional import Quadrature
from orpheus.numerics.symmetry import SubgroupOfO3
from tests.numerics import NUMERICS_DATA

ORPHEUS_ROOT = pathlib.Path(__file__).resolve().parents[2]
_NUMERICS = ORPHEUS_ROOT / "orpheus" / "numerics"

# ===========================================================================
# Helpers — every one `_r2`-prefixed, every rule built INSIDE a test body
# ===========================================================================

#: Labels only.  A production call in a parametrize list runs at import, i.e. AFTER a
#: mutation plugin installs, so a raising mutant kills COLLECTION and a `^FAILED`
#: scanner reads "0 caught" (vv Mode-8, third pipeline class).
_R2_RULE_LABELS: tuple[str, ...] = (
    "lebedev(5)", "lebedev(9)", "lebedev(13)", "lebedev(17)",
    "level_symmetric(4)", "level_symmetric(6)", "level_symmetric(8)",
    "product(4,8)", "product(8,8)",
    "gauss_legendre(2)", "gauss_legendre(8)", "gauss_legendre(16)",
    "product(4,4)", "folded_product(2,4)", "folded_product(4,8)",
)
_R2_FOLD_LABELS: tuple[str, ...] = (
    "folded_product(2,4)", "folded_product(4,4)", "folded_product(4,8)",
    "folded_product(6,8)", "folded_product(8,8)",
)


def _r2_rule(label: str) -> Quadrature:
    """Build a shipped rule from its LABEL — never from a parametrize list."""
    fam, _, args = label.partition("(")
    nums = [int(t) for t in args.rstrip(")").split(",") if t.strip()]
    if fam == "lebedev":
        return Quadrature.lebedev(nums[0])
    if fam == "level_symmetric":
        return Quadrature.level_symmetric(nums[0])
    if fam == "gauss_legendre":
        return Quadrature.gauss_legendre(nums[0])
    if fam == "product":
        return Quadrature.product(nums[0], nums[1])
    if fam == "folded_product":
        return Quadrature.folded_product(nums[0], nums[1])
    raise ValueError(f"unknown rule label {label!r}")


def _r2_groups() -> list[SubgroupOfO3]:
    """The 31 expressible members — the campaign's own denominator."""
    G = SubgroupOfO3
    out = [G.Trivial, G.Dinfh, G.OctahedralOh, G.IcosahedralIh, G.SO3, G.O3]
    out += [G.Mirror(a) for a in "xyz"]
    out += [G.SO2(a) for a in "xyz"]
    out += [G.O2(a) for a in "xyz"]
    out += [G.Cn(n) for n in range(1, 9)]
    out += [G.Dnh(n) for n in range(1, 9)]
    return out


def _r2_finite_groups() -> list[SubgroupOfO3]:
    return [g for g in _r2_groups() if g.realization.is_finite]


def _r2_continuous_groups() -> list[SubgroupOfO3]:
    return [g for g in _r2_groups() if not g.realization.is_finite]


def _r2_chart_spelling(measure: DiscreteMeasure) -> DiscreteMeasure:
    """The SAME orbit-space measure with its nodes stored at CHART width.

    ``Quotient.contains`` accepts both spellings, so one mathematical object has two
    storages — which is exactly the twin C3 names.
    """
    assert isinstance(measure.support, Quotient), "not an orbit-space measure"
    chart = np.asarray(
        measure.support.orbit_coordinates(measure.nodes), dtype=float
    )
    return DiscreteMeasure(
        nodes=chart, weights=measure.weights, support=measure.support
    )


def _r2_off_axis_sphere_measure(eps: float) -> DiscreteMeasure:
    """Two normalised nodes at :math:`(\\pm 1, \\varepsilon, 0)` on the BARE sphere.

    A bare support means the quotienting group is ``{e}`` — FINITE — so the kernel's
    position leg genuinely RUNS (on a ``S^2/O(2)_a`` support the barycentre lies on the
    axis by construction and the leg is a tautology).  The reference is arithmetic:
    the set is :math:`SO(2)_x`-closed iff :math:`\\varepsilon` is within the window.
    """
    pts = np.array([[1.0, eps, 0.0], [-1.0, eps, 0.0]])
    pts /= np.linalg.norm(pts, axis=1)[:, None]
    return DiscreteMeasure(
        nodes=pts, weights=np.array([2.0, 2.0]), support=SPHERE
    )


def _r2_two_azimuth_pairs() -> np.ndarray:
    r"""Four equatorial unit nodes at :math:`\varphi = 0,\;5\!\cdot\!10^{-10},\;
    \pi/2,\;\pi/2 + 5\!\cdot\!10^{-8}`.

    The first pair is below EVERY candidate window, the second sits between
    :math:`10^{-9}` (today's literal) and :math:`10^{-7}`.  So the count discriminates
    the literal from the argument, and nothing else does: `[M]` every shipped rule's
    azimuths are separated by at least :math:`2\pi/40`, i.e. :math:`10^8` above both
    windows.
    """
    phi = np.array([0.0, 5e-10, np.pi / 2.0, np.pi / 2.0 + 5e-8])
    return np.stack([np.cos(phi), np.sin(phi), np.zeros(4)], axis=1)


def _r2_asymmetric_sphere_measure() -> DiscreteMeasure:
    """The elegance review's own fixture — 9 seeded random unit nodes, distinct weights."""
    rng = np.random.default_rng(7)
    pts = rng.normal(size=(9, 3))
    pts /= np.linalg.norm(pts, axis=1)[:, None]
    return DiscreteMeasure(
        nodes=pts, weights=np.abs(rng.normal(size=9)) + 0.5, support=SPHERE
    )


class _R2OrderStub:
    """A pre-order element for :func:`_maximal` — the ONLY reachable witness.

    ``_maximal`` reads exactly two things off its items: ``!=`` and ``.contains``.  So a
    six-line duck type is its honest unit input, and it is the only way to build the pair
    the change is about: `[M]` over the 31 EXPRESSIBLE ``SubgroupOfO3`` members there are
    **0** distinct pairs that contain each other, so no production value can tell
    ``h != g`` from ``h`` strictly containing ``g`` (`vv` #14 / lessons L43d — split the
    pure predicate so a synthetic input is a two-line test).
    """

    def __init__(self, tag: str, below: "set[str]") -> None:
        self.tag = tag
        self._below = below

    def contains(self, other: "_R2OrderStub") -> bool:
        return other.tag in self._below or other.tag == self.tag

    def __eq__(self, other: object) -> bool:
        return isinstance(other, _R2OrderStub) and other.tag == self.tag

    def __hash__(self) -> int:
        return hash(self.tag)

    def __repr__(self) -> str:
        return f"_R2OrderStub({self.tag!r})"


def _r2_module_scope_imports(path: pathlib.Path) -> list[str]:
    """Every RUNTIME module-scope import target of ``path`` (TYPE_CHECKING excluded)."""
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    out: list[str] = []
    for node in tree.body:                       # module scope only
        if isinstance(node, ast.ImportFrom):
            base = "." * node.level + (node.module or "")
            out.append(base)
        elif isinstance(node, ast.Import):
            out.extend(a.name for a in node.names)
    return out


def _r2_deferred_imports(path: pathlib.Path) -> list[tuple[int, str]]:
    """Imports inside a function body (depth > 0) — ``(lineno, target)``."""
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    found: list[tuple[int, str]] = []

    class _V(ast.NodeVisitor):
        depth = 0

        def visit_FunctionDef(self, node):        # noqa: N802
            self.depth += 1
            self.generic_visit(node)
            self.depth -= 1

        def visit_AsyncFunctionDef(self, node):   # noqa: N802
            # spelled out rather than aliased: `visit_AsyncFunctionDef =
            # visit_FunctionDef` is a pyright signature error, and a
            # `# type: ignore` there would suppress the very static check that
            # catches this class of defect (coding-elegance anti-pattern #19).
            self.depth += 1
            self.generic_visit(node)
            self.depth -= 1

        def visit_Import(self, node):             # noqa: N802
            if self.depth:
                found.extend((node.lineno, a.name) for a in node.names)

        def visit_ImportFrom(self, node):         # noqa: N802
            if self.depth:
                found.append((node.lineno, "." * node.level + (node.module or "")))

    _V().visit(tree)
    return found


def _r2_call_sites(name: str, root: pathlib.Path) -> list[str]:
    """Every ``name(...)`` call site under ``root`` — ``file:lineno``."""
    out: list[str] = []
    for path in sorted(root.rglob("*.py")):
        if "__pycache__" in str(path):
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            fn = node.func
            spelled = (
                fn.id if isinstance(fn, ast.Name)
                else fn.attr if isinstance(fn, ast.Attribute)
                else None
            )
            if spelled == name:
                out.append(f"{path.relative_to(ORPHEUS_ROOT)}:{node.lineno}")
    return out


class _R2CountingSpy:
    """Wrap a module attribute, record every call, and ASSERT its own installation.

    A wrapping census that binds nothing prints a clean, confident zero — the
    safe-looking direction (agent memory, harness discipline).  So the constructor
    raises unless it rebound at least one module, and the count is read from the
    instance, never from a banner.
    """

    def __init__(self, monkeypatch, module, attr: str) -> None:
        self.calls: list[tuple] = []
        original = getattr(module, attr)
        bound = 0
        for candidate in list(sys.modules.values()):
            if getattr(candidate, attr, None) is original:
                monkeypatch.setattr(candidate, attr, self._wrap(original), raising=True)
                bound += 1
        if bound == 0:
            raise RuntimeError(
                f"the {attr!r} spy rebound 0 modules — the instrument is not "
                f"installed, and a zero it reports would be meaningless"
            )
        self.bound = bound

    def _wrap(self, original):
        def wrapper(*args, **kwargs):
            self.calls.append((args, kwargs))
            return original(*args, **kwargs)

        return wrapper


# ===========================================================================
# Group A — the import graph  (-> tests/test_layer_imports.py)
# ===========================================================================


class TestR2ImportGraph:
    r"""The import graph reads like the mathematics, and it is ENFORCED.

    ``geometry.transformation`` <- ``symmetry`` <- ``manifold`` <- ``measure`` <-
    ``invariance``.  The end-to-end half (A1) is a RE-KEY of the module's existing
    ``test_entry_point_imports_in_a_fresh_interpreter``: its parametrize list grows from
    6 entries to 10.

    ⛔ `[M]` injected and RUN on a shadow copy, one subprocess per (variant, entry)
    (``scratch/_r2_import_probe.py``): with ``manifold -> symmetry`` at module scope AND
    ``symmetry`` still importing ``AXIS_INDEX`` from ``manifold``, **6 of 9** entry points
    die with ``ImportError: cannot import name 'AXIS_INDEX' from partially initialized
    module``.  ``import orpheus`` alone still returns rc=0 — the order-dependence
    ``plan-authoring`` §6d warns about, which is why A2 exists beside A1: A1 is
    order-dependent, A2 is not.
    """

    @pytest.mark.foundation
    def test_a2_symmetry_does_not_import_manifold_or_measure_at_runtime(self) -> None:
        """The design-level gate: the 2-cycle is unspellable, not merely unreached.

        A1 can pass by import ORDER; this cannot.  ``TYPE_CHECKING`` references are
        allowed — they are erased at runtime — so the check reads module-scope
        ``Import``/``ImportFrom`` nodes only.
        """
        targets = _r2_module_scope_imports(_NUMERICS / "symmetry.py")
        forbidden = [
            t for t in targets
            if t.lstrip(".") in ("manifold", "measure",
                                 "orpheus.numerics.manifold",
                                 "orpheus.numerics.measure")
        ]
        assert not forbidden, (
            f"symmetry.py imports {forbidden} at module scope; with manifold's own "
            f"module-scope import of symmetry that is a direct 2-cycle "
            f"(measured: 6 of 9 entry points die)"
        )
        # positive control: the import census sees the edges that DO survive
        assert any("geometry.transformation" in t for t in targets), (
            "the AST census found no geometry.transformation import — the filter is "
            "broken, and its negative above carries no information"
        )

    @pytest.mark.foundation
    def test_a3_no_cycle_motivated_lazy_imports_remain(self) -> None:
        r"""No function-scope import of ``symmetry``/``manifold``/``measure``/
        ``invariance`` survives in the four modules.

        ⚠ The predicate is the CYCLE-motivated set, not "all deferred imports":
        `[M]` ``measure.py`` legitimately defers ``space``, ``axis``,
        ``generating_measure`` (×3) and ``exactness`` — six imports against modules R2
        does not restructure — and ``manifold.py`` defers ``sympy`` ×3 and ``LEGENDRE``.
        A done-when reading "only sympy/LEGENDRE/itertools" would be unreachable, and
        would be chased and then relaxed (``plan-authoring`` §10).
        """
        cycle_targets = {"symmetry", "manifold", "measure", "invariance"}
        offenders: list[str] = []
        for name in ("measure.py", "manifold.py", "symmetry.py", "invariance.py"):
            path = _NUMERICS / name
            for lineno, target in _r2_deferred_imports(path):
                tail = target.lstrip(".").split(".")[-1]
                if tail in cycle_targets:
                    offenders.append(f"{name}:{lineno} -> {target}")
        assert not offenders, (
            "cycle-motivated lazy imports must be gone after R2:\n  "
            + "\n  ".join(offenders)
        )

    @pytest.mark.foundation
    def test_a4_manifold_no_longer_duck_types_the_group(self) -> None:
        """``_trivial_group`` is retired and ``Quotient.by`` is a real annotation.

        ⚠ Paired with a POSITIVE control: a mis-typed module attribute also reads
        ``False`` under ``hasattr``, so the absence leg alone is unfalsifiable.
        """
        assert hasattr(man, "Quotient"), (
            "the positive control failed — this module is not manifold, so the "
            "absence assertions below carry no information"
        )
        assert not hasattr(man, "_trivial_group")
        # Resolve the field's annotation in manifold's OWN namespace: a name that
        # is imported under TYPE_CHECKING only is absent there and raises
        # (`typing.get_type_hints` would trip on the OTHER deferred names first —
        # `ReferenceMeasure` is still TYPE_CHECKING-only, legitimately).
        annotation = Quotient.__dataclass_fields__["by"].type
        resolved = eval(annotation, vars(man))
        assert resolved is SubgroupOfO3, (
            f"Quotient.by resolves to {resolved!r}, not the real class — the "
            f"TYPE_CHECKING forward reference is still in force"
        )


# ===========================================================================
# Group B — ONE closure  (-> tests/numerics/test_invariance.py)
# ===========================================================================


class TestR2OneClosure:
    r""""One closure" is a STRUCTURAL claim, so it is asserted structurally.

    Until R2 the claim was made by three docstrings and contradicted by the code:
    ``_invariance_on_orbit_space`` inlined a character-for-character copy of
    ``_orbit_space_closure``'s body (finding C1).  `[M]` at HEAD ``_orbit_closure`` has
    **2** call sites; after R2 it has **1**, so B1 is RED BEFORE the carve — which is
    what ``plan-authoring`` §6c asks of every gate.
    """

    @pytest.mark.foundation
    def test_b1_the_closure_has_exactly_one_caller(self) -> None:
        sites = [
            s for s in _r2_call_sites("_orbit_closure", ORPHEUS_ROOT / "orpheus")
        ]
        assert len(sites) == 1, (
            f"_orbit_closure is called from {len(sites)} sites {sites}; the "
            f"'one closure' claim is that there is exactly one"
        )

    @pytest.mark.foundation
    @pytest.mark.parametrize(
        "rule_label,group_name",
        [("product(4,8)", "Dnh8"), ("folded_product(4,8)", "sigma_y"),
         ("gauss_legendre(8)", "sigma_x")],
    )
    def test_b2_both_verbs_enter_the_same_closure_with_the_same_arguments(
        self, monkeypatch, rule_label: str, group_name: str
    ) -> None:
        """A COUNTING SPY, not a reading.

        Two verbs claiming to share a path is a claim about the ROUTE, and a route claim
        cannot be gated by an output (`vv` #26): a correct implementation and one that
        recomputes the match are indistinguishable in the answer.
        """
        measure = _r2_rule(rule_label).measure
        group = {
            "Dnh8": SubgroupOfO3.Dnh(8),
            "sigma_y": SubgroupOfO3.Mirror("y"),
            "sigma_x": SubgroupOfO3.Mirror("x"),
        }[group_name]

        spy = _R2CountingSpy(monkeypatch, inv, "_orbit_closure")
        measure.is_invariant_under(group)
        n_after_predicate = len(spy.calls)
        measure.certificate_under(group)
        n_after_certificate = len(spy.calls)

        assert n_after_predicate == 1, (
            f"is_invariant_under entered the closure {n_after_predicate} times"
        )
        assert n_after_certificate == 2, (
            f"certificate_under did not enter the closure "
            f"({n_after_certificate - n_after_predicate} entries)"
        )
        (a_args, _), (b_args, _) = spy.calls[0], spy.calls[1]
        np.testing.assert_array_equal(
            np.asarray(a_args[0]), np.asarray(b_args[0]),
            err_msg="the two verbs pass DIFFERENT node arrays to the one closure",
        )
        np.testing.assert_array_equal(
            np.asarray(a_args[1]), np.asarray(b_args[1]),
            err_msg="the two verbs pass DIFFERENT weight arrays",
        )

    @pytest.mark.foundation
    @pytest.mark.parametrize("rule_label", _R2_RULE_LABELS)
    def test_b3_predicate_equals_certificate_existence_for_finite_groups(
        self, rule_label: str
    ) -> None:
        """`[M]` 330 rows (22 finite groups x 15 shipped rules), 0 mismatches at HEAD."""
        measure = _r2_rule(rule_label).measure
        bad = [
            g.name for g in _r2_finite_groups()
            if measure.is_invariant_under(g)
            != (measure.certificate_under(g) is not None)
        ]
        assert not bad, f"{rule_label}: predicate and certificate disagree on {bad}"

    @pytest.mark.foundation
    def test_b3c_the_continuous_asymmetry_is_POLICY_and_is_pinned(self) -> None:
        r"""A continuous group has no finite element set, so it has no certificate —
        while it may perfectly well leave the measure invariant.

        `[M]` **135 rows** (9 continuous members x 15 rules): ``certificate_under`` is
        ``None`` on **135 of 135**, and ``is_invariant_under`` is ``True`` on exactly
        **6** — :math:`\{SO(2)_x, O(2)_x\}` on ``gauss_legendre(2/8/16)``.

        This row exists so a later "simplification" cannot make the certificate the
        kernel: that would silently report every polar marginal as NOT invariant under
        the very group its orbit space is named by.
        """
        witnesses: list[str] = []
        rows = 0
        for label in _R2_RULE_LABELS:
            measure = _r2_rule(label).measure
            for g in _r2_continuous_groups():
                rows += 1
                assert measure.certificate_under(g) is None, (
                    f"{label} x {g.name}: a continuous group produced a certificate"
                )
                if measure.is_invariant_under(g):
                    witnesses.append(f"{label} x {g.name}")
        assert rows == 135, rows
        assert sorted(witnesses) == sorted([
            "gauss_legendre(2) x SO2_x", "gauss_legendre(2) x O2_x",
            "gauss_legendre(8) x SO2_x", "gauss_legendre(8) x O2_x",
            "gauss_legendre(16) x SO2_x", "gauss_legendre(16) x O2_x",
        ]), witnesses

    @pytest.mark.foundation
    @pytest.mark.parametrize(
        "rule_label,group_name", [("product(4,8)", "sigma_z"),
                                  ("folded_product(4,8)", "sigma_y"),
                                  ("lebedev(9)", "sigma_x"),
                                  ("gauss_legendre(8)", "sigma_x")],
    )
    def test_b4_single_motion_face_agrees_with_the_group_face(
        self, rule_label: str, group_name: str
    ) -> None:
        """``permutation_under`` is the one-element reading of ``certificate_under``."""
        measure = _r2_rule(rule_label).measure
        group = SubgroupOfO3.Mirror(group_name[-1])
        cert = measure.certificate_under(group)
        assert cert is not None, f"{rule_label} is not {group.name}-invariant"
        for motion, perm in zip(cert.operators, cert.permutations):
            got = measure.permutation_under(motion.linear_part)
            assert got is not None
            np.testing.assert_array_equal(got.indices, perm.indices)


# ===========================================================================
# Group C — the deleted step 2  (-> tests/numerics/test_invariance.py)
# ===========================================================================


class TestR2DeletedStepTwo:
    r"""Kernel step 2 (:math:`G \subseteq H \Rightarrow` ``True``) was an OPTIMISATION,
    and the closure re-proves it.

    ⚠ Deleting it reddens nothing by itself (`[M]` qa's M5/M6 arms over 2290 kernel
    calls), which is exactly why the honest deliverable is this MUST-STAY-GREEN table
    rather than a mutation count.  The rows are ENUMERATED from the tree
    (:math:`H` = the support's quotienting group; every candidate it contains), never
    tabulated: `[M]` **32 distinct (rule x group) rows** — 10 bare-sphere rules x
    ``{Trivial}``, 3 ``gauss_legendre`` rules x
    ``{Trivial, sigma_y, sigma_z, SO2_x, O2_x, D_1h}``, 2 folds x ``{Trivial, sigma_y}``.

    ⛔ The kernel's docstring must SAY step 2 was an optimisation; otherwise the next
    reader restores it as a missing guard.
    """

    @pytest.mark.foundation
    @pytest.mark.parametrize("rule_label", _R2_RULE_LABELS)
    def test_c1_every_subgroup_of_the_stabiliser_is_still_invariant(
        self, rule_label: str
    ) -> None:
        measure = _r2_rule(rule_label).measure
        support = measure.support
        stabiliser = (
            support.by if isinstance(support, Quotient) else SubgroupOfO3.Trivial
        )
        rows = [
            g for g in _r2_groups()
            if g.normalises(stabiliser) and stabiliser.contains(g)
        ]
        assert rows, f"{rule_label}: the step-2 row set is empty — enumerate it again"
        failed = [g.name for g in rows if not measure.is_invariant_under(g)]
        assert not failed, (
            f"{rule_label}: the closure failed to re-prove step 2 for {failed}; "
            f"the deletion was NOT a theorem"
        )


# ===========================================================================
# Group D — the position test at the NODE window
# ===========================================================================


class TestR2PositionWindow:
    r"""The identity component's position test runs at the NODE window
    (``atol * _NODE_WINDOW_FACTOR``), not at the WEIGHT window.

    ⛔ **The change is INERT on every shipped rule** — `[M]` each rule's off-axis
    barycentre residual is either ``0.000e+00`` (bit-exact: the axial lift puts
    :math:`\mu` ON the axis) or ``>= 5.774e-01``, an :math:`10^{11}` separation, so
    **0 of 15** rules has a residual in :math:`(10^{-13}, 10^{-11}]`.  Its only witness
    is manufactured, and D2 records the inertness so a future audit does not read D1's
    green as production coverage.
    """

    @pytest.mark.foundation
    @pytest.mark.parametrize("eps,expected", [(1e-14, True), (1e-12, True), (1e-10, False)])
    def test_d1_admission_tracks_the_node_window(self, eps: float, expected: bool) -> None:
        """`[M]` at HEAD: 1e-14 True, 1e-12 **False**, 1e-10 False.  The middle row is
        the change's only witness; the outer two are its controls."""
        measure = _r2_off_axis_sphere_measure(eps)
        for axis_group in (SubgroupOfO3.SO2("x"), SubgroupOfO3.O2("x")):
            assert measure.is_invariant_under(axis_group, atol=1e-13) is expected, (
                f"eps={eps:.0e}: {axis_group.name} admission should be {expected} at "
                f"the node window 1e-11"
            )

    @pytest.mark.foundation
    def test_d2_the_window_change_is_inert_on_every_shipped_rule(self) -> None:
        """The separation table, asserted rather than believed."""
        in_band: list[str] = []
        for label in _R2_RULE_LABELS:
            bary = _embedded_nodes(_r2_rule(label).measure)
            for axis in range(3):
                off = float(np.abs(np.delete(bary, axis, axis=1)).max())
                if 1e-13 < off <= 1e-11:
                    in_band.append(f"{label}[axis {axis}] = {off:.3e}")
        assert not in_band, (
            "a shipped rule now sits in the widened band; the change is no longer "
            f"inert and D1's manufactured witness is no longer the only one: {in_band}"
        )

    @pytest.mark.foundation
    def test_d3_the_window_is_read_from_the_constant_not_respelled(
        self, monkeypatch
    ) -> None:
        """A literal ``1e-11`` would drift from ``_NODE_WINDOW_FACTOR`` silently."""
        seen: list[float] = []
        original = sym.IdentityComponent.fixes

        def spy(self, points, *, atol):
            seen.append(atol)
            return original(self, points, atol=atol)

        monkeypatch.setattr(sym.IdentityComponent, "fixes", spy, raising=True)
        measure = _r2_off_axis_sphere_measure(1e-12)
        measure.is_invariant_under(SubgroupOfO3.SO2("x"), atol=1e-13)
        assert seen, "the position leg never ran — the spy has nothing to report"
        expected = 1e-13 * inv._NODE_WINDOW_FACTOR
        assert all(a == pytest.approx(expected, rel=0, abs=0) for a in seen), (
            f"the position test used {seen}, not atol * _NODE_WINDOW_FACTOR "
            f"= {expected}"
        )


# ===========================================================================
# Group E — the azimuth window comes from the argument
# ===========================================================================


class TestR2AzimuthWindow:
    r"""``_distinct_azimuths`` takes its angular window from its ARGUMENT.

    Until R2 a literal ``1e-9`` shadowed the ``atol`` the function was handed (C4), so
    the function's own parameter could not move its answer.  ⛔ Inert on production:
    `[M]` ``n_az`` under merge ``1e-9`` and under merge ``1e-13`` is identical on
    **15 of 15** shipped rules.
    """

    @pytest.mark.foundation
    @pytest.mark.parametrize("atol,expected", [(1e-13, 4), (1e-9, 3), (1e-7, 2)])
    def test_e1_the_angular_window_is_the_argument(self, atol: float, expected: int) -> None:
        """`[M]` at HEAD **all three rows read 3** — the literal answers for every
        argument.  After R2 they read **4 / 3 / 2**, so TWO of the three are witnesses
        and the middle one is the control that pins the window the literal used to be.
        (The draft first predicted ``3`` at ``atol=1e-13``; the dry-run refuted it —
        with the merge window at ``1e-13`` the ``5e-10`` pair no longer merges either.)
        """
        got = _distinct_azimuths(_r2_two_azimuth_pairs(), atol)
        assert got == expected, (
            f"atol={atol:.0e}: expected {expected} distinct azimuths, got {got}; a "
            f"hard-coded merge window cannot see its own argument"
        )

    @pytest.mark.foundation
    def test_e2_the_tighter_window_moves_nothing_shipped(self) -> None:
        """The declared-blindness row, with the table it declares."""
        expected = {
            "lebedev(5)": 8, "lebedev(9)": 16, "lebedev(13)": 24, "lebedev(17)": 40,
            "level_symmetric(4)": 12, "level_symmetric(6)": 20, "level_symmetric(8)": 36,
            "product(4,8)": 8, "product(8,8)": 8,
            "gauss_legendre(2)": 2, "gauss_legendre(8)": 2, "gauss_legendre(16)": 2,
            "product(4,4)": 4, "folded_product(2,4)": 2, "folded_product(4,8)": 2,
        }
        got = {
            label: _distinct_azimuths(_embedded_nodes(_r2_rule(label).measure), 1e-13)
            for label in _R2_RULE_LABELS
        }
        assert got == expected, got


# ===========================================================================
# Group F — candidate_groups reads the EMBEDDED nodes
# ===========================================================================


class TestR2Candidates:
    r"""One measure, one candidate set — whatever width its nodes are stored at.

    Until R2 ``candidate_groups`` branched on ``measure.nodes.shape[1]``, so ONE fold had
    TWO answers by SPELLING (C3): `[M]` ``folded_product(4,8)`` stored at ambient width
    offered 20 candidates and reported :math:`\{D_{2h}\}`, the same object stored at
    chart width offered 15 and reported
    :math:`\{\sigma_x, \sigma_y, \sigma_z\}`.

    ⚠ Every assertion below compares SETS.  `[M]` the walk and the bruteforce scan return
    different tuple ORDER on a chart-spelled fold (``['sigma_z','sigma_x','sigma_y']`` vs
    ``['sigma_x','sigma_y','sigma_z']``), so a tuple comparison would be a flake.
    """

    @pytest.mark.foundation
    @pytest.mark.parametrize("fold_label", _R2_FOLD_LABELS)
    def test_f1_both_spellings_offer_the_same_candidates(self, fold_label: str) -> None:
        ambient = _r2_rule(fold_label).measure
        chart = _r2_chart_spelling(ambient)
        assert {g.name for g in candidate_groups(ambient)} == {
            g.name for g in candidate_groups(chart)
        }

    @pytest.mark.foundation
    @pytest.mark.parametrize("fold_label", _R2_FOLD_LABELS)
    def test_f2_because_the_barycentres_agree_bit_exactly(self, fold_label: str) -> None:
        """The premise F1 rests on — asserted, not assumed.  `[M]` ``max|diff|`` is
        ``0.000e+00`` on ``folded_product(4,8)``."""
        ambient = _r2_rule(fold_label).measure
        np.testing.assert_array_equal(
            _embedded_nodes(ambient), _embedded_nodes(_r2_chart_spelling(ambient))
        )

    @pytest.mark.foundation
    @pytest.mark.parametrize("fold_label", _R2_FOLD_LABELS)
    def test_f3_the_fold_reports_D2h_in_both_spellings_and_both_methods(
        self, fold_label: str
    ) -> None:
        ambient = _r2_rule(fold_label).measure
        for measure in (ambient, _r2_chart_spelling(ambient)):
            walk = {g.name for g in symmetry_groups(measure, method="walk")}
            brute = {g.name for g in symmetry_groups(measure, method="bruteforce")}
            assert walk == brute, (walk, brute)
            assert walk == {"D_2h"}, walk

    @pytest.mark.foundation
    @pytest.mark.parametrize("n", [2, 8, 16])
    def test_f4_a_polar_marginals_maximal_groups_are_O2x_and_D2h(self, n: int) -> None:
        r"""⚠ **RE-KEY of ``test_symmetry.py``'s slab row.** `[M]` at HEAD the answer is
        :math:`\{O(2)_x, \sigma_x\}`; after R2 the candidate set is derived from the same
        barycentres the invariance test reads, so :math:`D_{2h}` is offered, is
        invariant, and ABSORBS :math:`\sigma_x`.

        The derivation legs below are what make this a computed answer rather than a pin.
        """
        measure = Quadrature.gauss_legendre(n_ordinates=n).measure
        walk = {g.name for g in symmetry_groups(measure, method="walk")}
        brute = {g.name for g in symmetry_groups(measure, method="bruteforce")}
        assert walk == brute, (walk, brute)
        assert walk == {"O2_x", "D_2h"}, walk

        d2h, o2x = SubgroupOfO3.Dnh(2), SubgroupOfO3.O2("x")
        sigma_x = SubgroupOfO3.Mirror("x")
        assert d2h.contains(sigma_x), "D_2h must absorb sigma_x, else the answer is wrong"
        assert not o2x.contains(sigma_x), "sigma_x FLIPS the axis; O(2)_x cannot hold it"
        assert measure.is_invariant_under(sigma_x), "sigma_x is still invariant"
        assert not measure.is_invariant_under(SubgroupOfO3.Dnh(4)), (
            "D_4h's C_4^z does not normalise O(2)_x — it must be refused"
        )
        assert not measure.is_invariant_under(SubgroupOfO3.Dinfh), (
            "D_inf-h's rotations about z do not normalise O(2)_x"
        )

    @pytest.mark.foundation
    def test_f5_the_fold_candidate_names_are_the_eighteen(self) -> None:
        r"""⚠ **RE-KEY of ``test_the_candidate_SET_is_untouched_by_the_carve``.**
        `[M]` at HEAD the set is 20 names; the barycentre projection zeroes the ``y``
        column, so ``folded_product(4,8)``'s azimuth count falls 4 -> 2 and ``C_4`` /
        ``D_4h`` are no longer offered.  The literal below is written independently of
        ``candidate_groups``.
        """
        got = sorted(g.name for g in candidate_groups(Quadrature.folded_product(4, 8).measure))
        assert got == sorted([
            "C_2", "D_1h", "D_2h", "Dinfh", "Ih",
            "O2_x", "O2_y", "O2_z", "O3", "Oh",
            "SO2_x", "SO2_y", "SO2_z", "SO3", "Trivial",
            "sigma_x", "sigma_y", "sigma_z",
        ]), got


# ===========================================================================
# Group G — _maximal on STRICT containment
# ===========================================================================


class TestR2Maximal:
    r"""``_maximal`` drops an element only when another STRICTLY contains it.

    Its docstring has always said "not strictly contained in another member"; the body
    tested ``h != g and h.contains(g)``, which additionally drops BOTH members of an
    equal-but-distinct pair.  ⛔ `[M]` **0 of the 31 expressible members** has such a
    partner (R1's ``Cn(1) -> Trivial`` normalisation closed the only known one), so no
    production value can tell the two predicates apart — G2 declares that blindness and
    G1 supplies the only reachable witness.
    """

    @pytest.mark.foundation
    def test_g1_an_equal_but_distinct_pair_survives(self) -> None:
        a = _R2OrderStub("a", {"b"})
        b = _R2OrderStub("b", {"a"})
        assert a != b and a.contains(b) and b.contains(a)
        assert set(_maximal([a, b])) == {a, b}, (
            "an equal-but-distinct pair must BOTH survive: neither is STRICTLY "
            "contained in the other"
        )

    @pytest.mark.foundation
    def test_g1b_a_strict_pair_still_collapses(self) -> None:
        """The control: strict containment must still drop the smaller."""
        big = _R2OrderStub("big", {"small"})
        small = _R2OrderStub("small", set())
        assert set(_maximal([big, small])) == {big}

    @pytest.mark.foundation
    def test_g2_no_expressible_pair_can_tell_the_predicates_apart(self) -> None:
        """The declared blindness, measured — so G1's green is not read as production
        coverage."""
        groups = _r2_groups()
        pairs = [
            (a.name, b.name) for a in groups for b in groups
            if a != b and a.contains(b) and b.contains(a)
        ]
        assert not pairs, (
            f"an expressible equal-but-distinct pair now exists {pairs}; G1's "
            f"synthetic witness is no longer the only one and this docstring is stale"
        )

    @pytest.mark.foundation
    @pytest.mark.parametrize("method", ["walk", "bruteforce"])
    def test_g3_an_asymmetric_rule_reports_the_trivial_group(self, method: str) -> None:
        r"""⛔ **Already true at HEAD** — `[M]` both methods read ``(Trivial,)`` today.

        This row is a REGRESSION PIN on **R1's** ``Cn(1) -> Trivial`` normalisation, not
        a witness for R2's ``_maximal`` change: mutating ``_maximal`` back to ``!=``
        leaves it GREEN, because the invariant set is the single element ``{Trivial}``
        and a one-element list has no pair to disagree about.  Stated here so the
        battery's null on this row reads as EXPECTED rather than as a blind gate.
        """
        got = {g.name for g in symmetry_groups(_r2_asymmetric_sphere_measure(), method=method)}
        assert got == {"Trivial"}, got


# ===========================================================================
# Group H — the verbs moved; the old names are gone
# ===========================================================================


class TestR2RetiredNames:
    """Every name R2 deletes or moves is GONE from its old home.

    ⚠ Each absence leg is paired with a POSITIVE control importing a name that DOES
    survive from the SAME module — otherwise a typo in the module path produces the same
    ``ImportError`` and the gate is unfalsifiable.
    """

    @pytest.mark.foundation
    def test_h1_the_group_no_longer_answers_for_a_measure(self) -> None:
        assert hasattr(SubgroupOfO3, "contains"), "positive control failed"
        assert not hasattr(SubgroupOfO3, "is_invariant"), (
            "SubgroupOfO3.is_invariant survives as a facade; a facade keeps one "
            "deferred import alive and the boundary is not moved"
        )

    @pytest.mark.foundation
    @pytest.mark.parametrize(
        "name",
        ["is_invariant", "_check_invariance", "_invariance_on_orbit_space",
         "maximal_invariance_groups", "orbit_certificate", "singular_set",
         "induced_permutation",                      # the names they LEFT under
         "certificate_under", "singular_set_under", "permutation_under",
         "candidate_groups", "_embedded_nodes", "_orbit_closure",
         "_NODE_WINDOW_FACTOR"],
    )
    def test_h2_the_kernel_names_left_the_symmetry_module(self, name: str) -> None:
        assert not hasattr(sym, name), (
            f"symmetry.{name} still resolves; the kernel did not move"
        )

    @pytest.mark.foundation
    def test_h2_control_the_group_names_stayed(self) -> None:
        """The control for the eleven absence rows above."""
        for name in ("SubgroupOfO3", "Realization", "IdentityComponent",
                     "_INCOMMENSURATE_ANGLES", "_axis_vector"):
            assert hasattr(sym, name), f"symmetry.{name} vanished — it must not"


class TestR2MeasureVerbs:
    """The five verbs live on ``DiscreteMeasure`` and DELEGATE to ``invariance``."""

    @pytest.mark.foundation
    @pytest.mark.parametrize(
        "verb", ["is_invariant_under", "certificate_under", "permutation_under",
                 "singular_set_under", "symmetry_groups"],
    )
    def test_h3_the_measure_carries_the_verb(self, verb: str) -> None:
        assert hasattr(DiscreteMeasure, verb)

    @pytest.mark.foundation
    def test_h3b_the_verb_delegates_rather_than_reimplements(self, monkeypatch) -> None:
        spy = _R2CountingSpy(monkeypatch, inv, "is_invariant_under")
        measure = Quadrature.product(4, 8).measure
        measure.is_invariant_under(SubgroupOfO3.Dnh(8))
        assert len(spy.calls) == 1, (
            f"the method entered invariance.is_invariant_under {len(spy.calls)} "
            f"times; a second implementation is a Pattern-2 twin"
        )

    @pytest.mark.foundation
    def test_h3c_quotient_reads_the_measures_own_certificate(self, monkeypatch) -> None:
        """``DiscreteMeasure.quotient`` goes through ``self.certificate_under``, not
        through a second (function-scope) import of the free function.

        ⚠ The spy targets the METHOD, not ``invariance.certificate_under``: the claim is
        about which surface ``quotient`` reads, and a spy on the free function would be
        satisfied by the retired lazy import too.
        """
        seen: list[object] = []
        original = DiscreteMeasure.certificate_under

        def spy(self, group, **kw):
            seen.append(group)
            return original(self, group, **kw)

        monkeypatch.setattr(DiscreteMeasure, "certificate_under", spy, raising=True)
        Quadrature.product(4, 8).measure.quotient(SubgroupOfO3.Mirror("y"))
        assert seen, (
            "quotient() did not read self.certificate_under — it still resolves the "
            "certificate through its own import"
        )


class TestR2OrdinatePermutationIsUnchanged:
    """H4b — the rename is a PURE re-spelling, so the answer is bit-identical.

    -> ``tests/numerics/test_quadrature_directional.py``.  The reference is a table
    FROZEN at HEAD before the carve (``tests/numerics/data/r2_perm_baseline.npz``, `[M]` 45
    (rule x axis) rows, 0 of them ``None``) — a value R2 cannot have produced.  Any
    tolerance here would be a red flag (`vv` bit-identity).
    """

    @pytest.mark.foundation
    @pytest.mark.parametrize("rule_label", _R2_RULE_LABELS)
    def test_h4b_every_permutation_is_array_equal_to_the_frozen_table(
        self, rule_label: str
    ) -> None:
        baseline = np.load(NUMERICS_DATA / "r2_perm_baseline.npz")
        quad = _r2_rule(rule_label)
        for axis in ("x", "y", "z"):
            motion = SelfPairedDeck.mirror(axis=axis, dimension=3).motion
            pi = quad.ordinate_permutation(motion)
            want = baseline[f"{rule_label}|{axis}"]
            if want.size == 1 and want[0] == -1:
                assert pi is None
            else:
                assert pi is not None
                np.testing.assert_array_equal(np.asarray(pi.indices), want)


class TestR2SelectionIsUnchanged:
    """H4a — quadrature SELECTION is the only production consumer of the predicate.

    -> ``tests/numerics/test_registry.py``.  ``AngularSymmetry.admits_symmetry``
    (``registry.py:1012``) is the ONE production call site of ``is_invariant``; R2
    re-spells it ``measure.is_invariant_under(self.owed)``.  A pure re-spelling, so the
    verdict must be EQUAL — chosen spec, its parameters, the point count, and every
    rejection STRING (the strings are an API the moment a test pins one) —
    except stage 0's, pinned at its STAGE since R3 of #434 re-worded it; its
    wording is pinned once in ``tests/numerics/test_registry.py`` (see
    ``stage_of``).

    The reference is a table FROZEN at HEAD before the carve
    (``tests/numerics/data/r2_selection_baseline.json``, `[M]` 48 (geometry x degree)
    verdicts, all
    48 selecting).
    """

    @pytest.mark.foundation
    @pytest.mark.parametrize(
        "geometry", ["slab", "sphere", "cylinder", "cartesian2d"]
    )
    def test_h4a_every_verdict_is_equal_to_the_frozen_table(self, geometry: str) -> None:
        import json

        from orpheus.numerics.quadrature.registry import select_quadrature

        baseline = json.loads((NUMERICS_DATA / "r2_selection_baseline.json").read_text())

        def stage_of(rejection: str) -> str:
            # Every rejection string is pinned VERBATIM except stage 0's,
            # which R3 of #434 re-worded (`[M]` 96 of 96 domain-mismatch
            # strings moved, 8 of 8 symmetry-mismatch strings unchanged, 0 of
            # 48 choices moved): that one is pinned at its STAGE here and its
            # wording is pinned once, deliberately, beside the R3 rows in
            # tests/numerics/test_registry.py.
            spec, reason = rejection.split("::", 1)
            if not reason.startswith("domain mismatch"):
                return rejection
            return f"{spec}::domain mismatch"

        for degree in (1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 17):
            key = f"{geometry}|{degree}"
            want = dict(baseline[key])
            want["rejected"] = sorted(stage_of(r) for r in want["rejected"])
            measure, log = select_quadrature(geometry, degree)
            got = {
                "chosen": log.chosen_spec.name if log.chosen_spec else None,
                "parameters": dict(log.chosen_parameters or {}),
                "n_points": int(measure.n_points),
                "rejected": sorted(stage_of(f"{n}::{r}") for n, r in log.rejected),
            }
            assert got == want, f"{key}: selection moved\n  want {want}\n  got  {got}"


# ===========================================================================
# Group I — C6: the ERR-045 message on a fold
# ===========================================================================


class TestR2Err045NamesTheOrbitSpace:
    r"""On a :math:`\sigma_y`-folded rule the specular guard must name the ORBIT SPACE.

    ⛔ **Not tagged ``catches("ERR-045")``.**  The guard fires today too — `[M]` on
    ``folded_product(4, 8)``, axis ``'y'``, the permutation IS the identity
    (``array_equal(perm, arange(16))``) because :math:`\sigma_y` acts trivially on
    :math:`S^2/\sigma_y`, and the raise reads *"same sign class instead of the outflow
    partner (ERR-045: wrong pairing, or a non-axis-aligned reflection that needs a
    different BC type)"*.  Neither diagnosis is the truth: the fold SPENT
    :math:`\sigma_y`, so the outflow half is not stored at all.  R2 changes the
    DIAGNOSIS, not the refusal, so a ``catches`` marker here would be a false coverage
    claim and a bare ``pytest.raises`` would be teeth-less — **only the two message legs
    attribute the change**.
    """

    @pytest.mark.foundation
    def test_i1_the_fold_refusal_names_the_spent_symmetry(self) -> None:
        quad = Quadrature.folded_product(4, 8)
        with pytest.raises(ReflectionDidNotMapInflowToOutflowError) as excinfo:
            assert_specular_pairing_maps_inflow_to_outflow(
                quad, "y", law_key="r2-probe"
            )
        message = str(excinfo.value)
        assert "SPENT the mirror" in message, (
            f"the refusal must name the fold that consumed the mirror; got: {message}"
        )
        assert "non-axis-aligned" not in message, (
            "the old misdiagnosis survives: the pairing is not wrong and the "
            "reflection is axis-aligned — the rule simply no longer stores the "
            "outflow half"
        )

    @pytest.mark.foundation
    def test_i2_control_the_unfolded_sibling_does_not_raise(self) -> None:
        """`[M]` no raise at HEAD — the guard is correct on an unfolded rule."""
        assert_specular_pairing_maps_inflow_to_outflow(
            Quadrature.product(4, 8), "y", law_key="r2-probe"
        )

    @pytest.mark.foundation
    def test_i3_control_the_folded_rules_other_axis_does_not_raise(self) -> None:
        """`[M]` no raise at HEAD — the ``x`` partner survives the ``y`` fold."""
        assert_specular_pairing_maps_inflow_to_outflow(
            Quadrature.folded_product(4, 8), "x", law_key="r2-probe"
        )
