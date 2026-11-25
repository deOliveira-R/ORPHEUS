r"""G03q: the stated face of a dense form is content (N03, spec 7.4; the user's ruling of 2026-11-24).

On the 21 harness Grams (17 frame-built, 4 hand-written), each built once in this process as the frame's own spelling,
``d = DenseMetric(DenseMetric.symmetrised(G), kernel=k, face=Face.GRAM)``: ``dataclasses.replace(d)`` keeps both faces' bytes; ``d.flipped().flipped()`` has the
stated matrix and the face of ``d``; ``d.flipped()``'s ``matrix`` and ``inverse_matrix`` are ``d``'s ``inverse_matrix`` and ``matrix``; the form stated as
``G+`` (``DenseMetric(d.matrix, kernel=k, face=METRIC)``) has a content digest different from ``d``'s (21 of 21: a field differs) and the same ``metric``
face bytes; one stated matrix at the two faces is two digests; ``replace(d, inverse_matrix=)`` and ``replace(d, matrix=)`` raise (both are derived);
``inverse_of`` is not an attribute of ``DenseMetric`` and ``DenseMetric(m, kernel=k)`` without ``face`` raises ``TypeError``; a matrix stated at either
face is taken as stored (``.matrix`` has the bytes of the matrix handed, the other face is its exact-count inverse); and, by ``ast``, no ``replace(`` in
``orpheus/`` passes a dense form's own field (``matrix``, ``inverse_matrix``, ``stated``, ``face`` or ``kernel``). The ``ast`` leg reads keywords, not
the receiver's type, so a ``replace`` on a form that passes none of the five fields is outside it.

Exact on any host: every byte comparison is of a computation with the same computation (``replace``, a double flip) or a field with itself; the gate
does not assert the host-arithmetic count of Grams on which ``replace`` moved the gram face before the unit. First red at ``ff305e0f``: import-only
(``Face``, ``flipped`` and the ``face`` keyword do not exist). Mutations: the second face re-derived from ``matrix`` (the N01 defect: the flip is not the other
face), ``flipped()`` returns ``self``, ``face`` dropped from the content, ``inverse_of`` kept as a method.
"""

from __future__ import annotations

import ast
import dataclasses
from pathlib import Path

import numpy as np
import pytest

from tests.gates.cross_method.w1 import _fixtures as fx

pytestmark = pytest.mark.foundation

_FIELDS = {"matrix", "inverse_matrix", "stated", "face", "kernel"}


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _grams() -> list[tuple[str, np.ndarray, object]]:
    out = []
    for family, _text, build in fx.frame_families():
        gram = np.asarray(build(), dtype=float)
        kernel, _ = fx.declared_kernel_of(family, gram)
        out.append((family, gram, kernel))
    return out


def test_the_stated_face_is_content() -> None:
    from orpheus.numerics.declaration import Coord
    from orpheus.numerics.metric import DenseMetric, Face
    from orpheus.values.content import content_digest

    grams = _grams()
    _require(len(grams) == 21, f"{len(grams)} Grams, expected the 21 of the harness")
    for family, gram, k in grams:
        d = DenseMetric(DenseMetric.symmetrised(gram), kernel=k, face=Face.GRAM)
        r = dataclasses.replace(d)
        _require(d.matrix.tobytes() == r.matrix.tobytes() and d.inverse_matrix.tobytes() == r.inverse_matrix.tobytes(), f"[{family}] replace moved a face")
        ff = d.flipped().flipped()
        _require(ff.stated.tobytes() == d.stated.tobytes() and ff.face is d.face, f"[{family}] two flips are not the identity")
        f = d.flipped()
        _require(f.face is not d.face, f"[{family}] flipped() kept the face")
        _require(f.matrix.tobytes() == d.inverse_matrix.tobytes() and f.inverse_matrix.tobytes() == d.matrix.tobytes(), f"[{family}] flipped() is not the other face")
        as_metric = DenseMetric(d.matrix, kernel=k, face=Face.METRIC)
        _require(content_digest(as_metric) != content_digest(d), f"[{family}] the form stated as G+ has the digest of inverse_of(G)'s form")
        _require(as_metric.matrix.tobytes() == d.matrix.tobytes(), f"[{family}] the metric face of the form stated as G+ differs")
        _require(content_digest(DenseMetric(d.stated, kernel=k, face=Face.METRIC)) != content_digest(DenseMetric(d.stated, kernel=k, face=Face.GRAM)), f"[{family}] the face is not content")
        # a matrix stated at either face is taken as stored; the other face is the exact-count inverse
        m = np.ascontiguousarray(d.stated)
        at_metric = DenseMetric(m, kernel=k, face=Face.METRIC)
        at_gram = DenseMetric(m, kernel=k, face=Face.GRAM)
        _require(at_metric.matrix.tobytes() == m.tobytes() and at_gram.inverse_matrix.tobytes() == m.tobytes(), f"[{family}] a stated matrix was not taken as stored")
        _require(at_metric.inverse_matrix.tobytes() == at_gram.matrix.tobytes(), f"[{family}] the second faces of the two statements differ")
    _, _, k0 = grams[0]
    d0 = DenseMetric(DenseMetric.symmetrised(grams[0][1]), kernel=k0, face=Face.GRAM)
    for kw in ("inverse_matrix", "matrix"):
        try:
            dataclasses.replace(d0, **{kw: d0.matrix})
        except (TypeError, ValueError):
            pass
        else:
            pytest.fail(f"replace(d, {kw}=) was accepted: the field is derived")
    _require(not hasattr(DenseMetric, "inverse_of"), "DenseMetric still has inverse_of")
    try:
        DenseMetric(np.eye(3), kernel=Coord(3, ()))  # type: ignore[call-arg]
    except TypeError:
        pass
    else:
        pytest.fail("DenseMetric(m, kernel=k) without a face was accepted")


def test_no_replace_in_orpheus_passes_a_dense_forms_field() -> None:
    import orpheus

    root = Path(orpheus.__file__).parent
    files = sorted(root.rglob("*.py"))
    _require(len(files) > 200, f"the census reads {len(files)} files")
    hits = []
    n_replace = 0
    for p in files:
        for n in ast.walk(ast.parse(p.read_text())):
            if isinstance(n, ast.Call) and ((isinstance(n.func, ast.Name) and n.func.id == "replace") or (isinstance(n.func, ast.Attribute) and n.func.attr == "replace" and isinstance(n.func.value, ast.Name) and n.func.value.id == "dataclasses")):
                n_replace += 1
                if any(k.arg in _FIELDS for k in n.keywords):
                    hits.append(f"{p.relative_to(root)}:{n.lineno}")
    _require(n_replace > 20, f"the census saw {n_replace} replace calls: the filter dropped them")
    _require(not hits, f"{len(hits)} of {n_replace} replace calls pass a dense form's field: {hits}")
