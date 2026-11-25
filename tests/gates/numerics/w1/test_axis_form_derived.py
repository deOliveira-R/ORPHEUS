r"""G03d, G03h, G03i: the axis carries one form from one home, the origin that minted it (N03, spec 7.4).

G03d. A weight and a form are not arguments of an axis on any of the five arms (``FromMeasure``, ``FromQuadrature``, ``FromBasis``, ``FromFrame``,
``Declared``): ``Axis(label, shape, weights=w, origin=...)`` raises ``TypeError`` for ``w`` equal to what the origin determines and for ``w``
disagreeing, a third positional argument raises, and ``Axis(..., form=DiagonalMetric(w), origin=...)`` raises ``TypeError``; ``Axis.of_measure`` and
``Axis.of_quadrature`` carry the measure's weights in ``axis.form``; the retired guard tag does not appear in ``numerics/axis.py``, the wave-tagged
issues of the ledger gate do not name it, and no ``determined_weights`` definition remains in ``orpheus/``.

G03h. ``Axis.restricted((0, 2))`` has shape ``(2,)``, the principal submatrix as its form, a ``restriction`` record (the parent's content digest and the
indices), origin ``Declared(NODAL, form)``, and refusals on a rank-2 axis, a ``Modal`` parent, an unsorted set, a repeated index and an empty set (fragment
"sorted ascending and unique") and on a dense parent (a declared edge tagged ``SCOPE-BOUNDARY[guard]``); a restricted axis and a ``Declared`` axis of equal
content compare unequal.

G03i. ``Axis.of_quadrature(q, "angular", face_axis=k)`` has form ``Diag(abs(cosine) * w + 0.0)`` for ``k = 0, 1, 2`` with zeros exactly where the cosine is
``0.0``, an identity different from the bulk axis, ``face_axis=3`` refused with the accessor's own error, and on a 1-D rule ``face_axis=1`` is an all-kernel
form by design.

Exact on any host (each comparison is of a computation with the computation it must equal; the counts are exact). First red at ``ff305e0f``: G03d is
behavioural (a ``weights=`` equal to the measure's is accepted: red against ``TypeError``; a disagreeing one is a ``ValueError``; the ``form=`` call is
already a ``TypeError``: green; the retired tag appears once); G03h is import-only on ``Axis.restricted``; G03i's ``face_axis=`` is a ``TypeError``.
Mutations: re-add the ``weights`` and ``form`` keywords, accepted and ignored (G03d); drop ``restriction`` from the identity (G03h); the cosine in place of
``abs(cosine)`` (G03i: the gate's own ``array_equal`` reddens, no refusal intervenes).
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import numpy as np
import pytest

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _orpheus_root() -> Path:
    import orpheus

    return Path(orpheus.__file__).parent


def _arms() -> list[tuple[str, tuple[str, tuple[int, ...]], object, np.ndarray]]:
    """``(arm, (label, shape), origin, weights the origin determines or ones)`` for the five arms."""
    from orpheus.numerics.axis import BasisKind, Declared, FromBasis, FromFrame, FromMeasure, FromQuadrature
    from orpheus.numerics.basis.spherical_harmonic_basis import SphericalHarmonicBasis
    from orpheus.numerics.frame_faces import angular_frame, basis_space
    from orpheus.numerics.quadrature.directional import Quadrature

    q = Quadrature.gauss_legendre(4)
    m = q.measure
    b = SphericalHarmonicBasis(L=1)
    fr = angular_frame(Quadrature.level_symmetric(4), 2)
    head = basis_space(fr).axes[0]
    shape_b = (int(b.layout.shape[0]),)
    return [
        ("FromMeasure", ("x", (4,)), FromMeasure(m), np.asarray(m.weights, dtype=float)),
        ("FromQuadrature", ("angular", (4,)), FromQuadrature(q), np.asarray(m.weights, dtype=float)),
        ("FromBasis", (b.layout.label, shape_b), FromBasis(b), np.ones(shape_b)),
        ("FromFrame", (head.label, head.shape), FromFrame(fr), np.ones(head.shape)),
        ("Declared", ("s", (3,)), Declared(BasisKind.NODAL), np.ones(3)),
    ]


def test_no_weights_or_form_argument_on_any_arm() -> None:
    from orpheus.numerics.axis import Axis
    from orpheus.numerics.metric import DiagonalMetric

    arms = _arms()
    _require(len(arms) == 5, f"{len(arms)} arms built, expected the five")
    for arm, (label, shape), origin, w in arms:
        for name, kw in (("equal weights", {"weights": w}), ("disagreeing weights", {"weights": 2.0 * w + 1.0}), ("form=", {"form": DiagonalMetric(w)})):
            try:
                Axis(label, shape, origin=origin, **kw)  # type: ignore[arg-type]
            except TypeError:
                continue
            except Exception as e:
                pytest.fail(f"[{arm}] {name}: raised {type(e).__name__} ({str(e)[:80]}), not TypeError")
            pytest.fail(f"[{arm}] {name}: the axis was constructed")
        try:
            Axis(label, shape, w, origin=origin)  # type: ignore[misc]
        except TypeError:
            pass
        else:
            pytest.fail(f"[{arm}] a third positional argument was accepted")


def test_the_generator_arms_carry_their_measure_in_the_form() -> None:
    from orpheus.numerics.axis import Axis
    from orpheus.numerics.metric import DiagonalMetric
    from orpheus.numerics.quadrature.directional import Quadrature

    q = Quadrature.gauss_legendre(4)
    w = np.asarray(q.measure.weights) + 0.0
    for name, ax in (("of_measure", Axis.of_measure(q.measure, "x")), ("of_quadrature", Axis.of_quadrature(q))):
        _require(isinstance(ax.form, DiagonalMetric), f"{name}: the form is {type(ax.form).__name__}")
        _require(np.asarray(ax.form.weights).tobytes() == np.ascontiguousarray(w).tobytes(), f"{name}: the form's weights are not the measure's")


def test_the_retired_guard_and_its_properties_are_gone() -> None:
    from tests.gates import test_elegance_debt_is_tagged as ledger

    axis_src = (_orpheus_root() / "numerics" / "axis.py").read_text()
    _require("ELEGANCE-DEBT[guard] #710" not in axis_src, "numerics/axis.py still carries the retired guard's tag")
    _require(710 not in ledger._WAVE_TAGGED_ISSUES, "the ledger gate's wave-tagged issues still name the retired guard")
    files = sorted(_orpheus_root().rglob("*.py"))
    _require(len(files) > 200, f"the census reads {len(files)} files")
    defs = [
        f"{p.relative_to(_orpheus_root())}:{n.lineno}"
        for p in files
        for n in ast.walk(ast.parse(p.read_text()))
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)) and n.name == "determined_weights"
    ]
    _require(not defs, f"{len(defs)} definitions of determined_weights remain: {defs}")
    _require(re.search(r"\bdetermined_weights\b", axis_src) is None, "numerics/axis.py still names determined_weights")


def test_axis_restricted() -> None:
    from orpheus.numerics.axis import Axis, BasisKind, Declared, Modal
    from orpheus.numerics.metric import DenseMetric, DiagonalMetric, Face
    from orpheus.numerics.declaration import Coord
    from orpheus.numerics.quadrature.directional import Quadrature
    ax = Axis.of_quadrature(Quadrature.gauss_legendre(4))
    r = ax.restricted((0, 2))
    _require(r.shape == (2,), f"the restricted axis has shape {r.shape}")
    _require(np.array_equal(r.form.weights, np.asarray(ax.form.weights)[[0, 2]]), "the form is not the principal submatrix")
    _require(r.restriction is not None and r.restriction.parent_key == ax.content_digest and tuple(r.restriction.indices) == (0, 2), f"the restriction record is {r.restriction}")
    _require(isinstance(r.origin, Declared) and r.origin.kind is BasisKind.NODAL, f"the origin is {r.origin}")
    refusals = (
        ("unsorted", lambda: ax.restricted((2, 0))),
        ("repeated", lambda: ax.restricted((1, 1))),
        ("empty", lambda: ax.restricted(())),
        ("rank 2", lambda: Axis("r2", (2, 2), origin=Declared(BasisKind.NODAL)).restricted((0,))),
        ("modal", lambda: Axis("m", (3,), origin=Declared(Modal(np.ones(3)))).restricted((0,))),
        (
            "dense parent",
            lambda: Axis("d", (3,), origin=Declared(BasisKind.NODAL, DenseMetric(np.diag([1.0, 2.0, 3.0]), kernel=Coord(3, ()), face=Face.METRIC))).restricted((0, 1)),
        ),
    )
    for label, call in refusals:
        try:
            call()
        except ValueError as e:
            if label in ("unsorted", "repeated"):
                _require("sorted ascending and unique" in str(e), f"[{label}] the refusal does not say 'sorted ascending and unique': {e}")
        else:
            pytest.fail(f"[{label}] restricted accepted it")
    declared = Axis("angular", (2,), origin=Declared(BasisKind.NODAL, DiagonalMetric(np.asarray(ax.form.weights)[[0, 2]])))
    _require(r != declared and r.content_digest != declared.content_digest, "a restricted axis equals a Declared axis of equal content")


def test_the_dense_parent_edge_is_tagged() -> None:
    from tests.gates import test_elegance_debt_is_tagged as ledger

    src = (_orpheus_root() / "numerics" / "axis.py").read_text()
    lines = src.splitlines()
    sites = [i for i, ln in enumerate(lines) if "SCOPE-BOUNDARY[guard]" in ln]
    _require(bool(sites), "numerics/axis.py carries no SCOPE-BOUNDARY[guard] (the dense-parent edge of Axis.restricted)")
    _require(not ledger._boundary_violations_in(src, label="axis.py"), "a boundary tag of axis.py lacks machinery, ruling or revisit")


def test_face_weighted_quadrature_axis() -> None:
    from orpheus.numerics.axis import Axis
    from orpheus.numerics.quadrature.directional import Quadrature

    q = Quadrature.product(2, 4)
    w = np.asarray(q.measure.weights)
    for k in range(3):
        c = np.asarray(q.mean_axis_cosine(k))
        a = Axis.of_quadrature(q, "angular", face_axis=k)
        fw = np.asarray(a.form.weights)
        _require(np.array_equal(fw, np.abs(c) * w + 0.0), f"face_axis={k}: the form is not abs(cosine) * w")
        _require(bool(np.array_equal(c == 0.0, fw == 0.0)), f"face_axis={k}: the zeros of the form are not the zeros of the cosine")
    _require(Axis.of_quadrature(q, "angular", face_axis=0) != Axis.of_quadrature(q), "the face-weighted axis has the identity of the bulk axis")
    try:
        Axis.of_quadrature(q, "angular", face_axis=3)
    except ValueError:
        pass
    else:
        pytest.fail("face_axis=3 was accepted")
    z = Axis.of_quadrature(Quadrature.gauss_legendre(4), "angular", face_axis=1)
    _require(bool((np.asarray(z.form.weights) == 0.0).all()), "a suppressed axis (face_axis=1 on a 1-D rule) is not an all-kernel form")
