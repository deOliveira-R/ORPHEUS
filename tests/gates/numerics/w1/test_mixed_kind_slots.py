r"""G04d: the kind of a slot is per slot, so the operator-data spaces of the transport equation construct (N04, spec 7.5).

The fission spectrum space ``chi = <energy_out FUNCTION>``, the scattering matrix ``<energy_out FUNCTION, energy_in MEASURE>`` and the fission dyad
``<energy_out FUNCTION, energy_in MEASURE, spatial (MEASURE, DENSITY)>`` construct (the operator-data types are ``[R]`` in the spec; the energy axes are formless, so their
variance is ``FormlessVariance``); ``dual()`` flips every slot's kind independently; ``dual(on={"energy_in"})`` moves that slot only.

First red at the base tip: import-only (``TensorSpace`` is absent). Mutation: one kind bit for the whole space as written raises at construction (an exception, so the
construction sits inside the gate's body and fails it); the replacement that reddens by assertion: ``dual()`` flips every slot to the first slot's flipped kind, reddens
"dual does not flip per slot". Exact on any host.
"""

from __future__ import annotations

import pytest

from tests.gates.numerics.w1 import _tensor_probe as P

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _spaces() -> dict[str, object]:
    from orpheus.numerics import factor as F

    TS = P.tensor_space()
    out_f = P.formless_factor(P.formless_axis(2, "energy_out"), F.FormlessVariance.FUNCTION)
    in_m = P.formless_factor(P.formless_axis(2, "energy_in"), F.FormlessVariance.MEASURE)
    sp = P.density_factor(P.diag_axis(P.SPATIAL_W, "spatial"))
    return {"chi": TS.of_axes(out_f), "scattering": TS.of_axes(out_f, in_m), "dyad": TS.of_axes(out_f, in_m, sp)}


def test_mixed_kind_spaces_construct_and_dual_flips_per_slot() -> None:
    for name, V in _spaces().items():
        D = V.dual()
        for f, g in zip(V.factors, D.factors):
            _require(g == f.dual(), f"[{name}] dual does not flip per slot: {f.axis.label} {f.variance} -> {g.variance}")
        kinds = {f.kind for f in V.factors}
        if len(kinds) > 1:
            _require(len({g.kind for g in D.factors}) > 1, f"[{name}] dual does not flip per slot: the mixed kinds collapsed to one")


def test_dual_on_moves_one_slot() -> None:
    V = _spaces()["dyad"]
    moved = V.dual(on={"energy_in"})
    _require(moved.factor("energy_in") == V.factor("energy_in").dual(), "dual(on=energy_in) did not flip energy_in")
    for l in ("energy_out", "spatial"):
        _require(moved.factor(l) == V.factor(l), f"dual(on=energy_in) moved {l}")
