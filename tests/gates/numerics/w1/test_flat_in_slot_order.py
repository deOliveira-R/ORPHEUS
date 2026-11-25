r"""G04e and G04f: the Kronecker walk is in slot order, and the energy window reads today's norm (N04, spec 7.5; law L8; decision ``[ON04-5]``).

G04e. On 200 seeded draws the gate spells the per-slot products itself (``_tensor_probe.slot_order_flat``: one broadcast product per formed slot, in slot order, or reversed, from the weights the FIXTURE stated and not from production's form objects), and
``flat_on(x, formed labels)`` is byte-equal to the slot-order spelling on 200 of 200 draws while the reversed spelling differs from it on at least 190 of the 200 (the opening measured
200 of 200 on the SN slab angular bulk space; 197 of 200 on the signed pair; the bound 190 is the gate's positive control that the instrument sees order), over (a) a positive angular-energy-spatial
space and (b) a signed pair ``Diag(1.3, -0.7, 2.1)`` and ``Diag(0.3, -1.1, 0.9, 0)`` (the second taken modulo its kernel). The harness' stored bytes (G00a) are the cross-check.
G04f. On a space ``<energy NoForm, spatial>`` and a seeded ``q``: ``norm(q)`` is bytewise ``sqrt(sum(q * apply_metric(q)))``; ``euclidean().norm(q)`` equals it; ``counting().norm(q)``
is bytewise ``sqrt(sum(q * q))`` (not ``np.linalg.norm``, a different last bit) and differs from ``norm(q)``.

First red at the base tip: import-only (``TensorSpace`` and ``flat_on`` absent; ``euclidean`` and ``counting`` absent). Mutations (described): reversing the walk in ``flat_on`` reddens the equality leg
of G04e (the gate's own slot-order spelling then differs from ``flat_on``; ``[R]`` the build runs it on the real ``TensorSpace``); ``euclidean`` returning ``counting()`` reddens G04f (norm 3.6237 against 4.0433 on the opening's draw);
``counting`` keeping the angular form reddens on a non-uniform vector. Exact on any host: the arithmetic is elementwise products and one sum in one process.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.gates.numerics.w1 import _tensor_probe as P

pytestmark = pytest.mark.foundation

N_DRAWS = 200
MIN_REVERSED_DIFFERS = 190


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _walk_rows(space: object, labels: tuple[str, ...], stated: dict[str, tuple[float, ...]]) -> tuple[int, int]:
    rng = np.random.default_rng(20261125)
    equal_slot = differs_reversed = 0
    for _ in range(N_DRAWS):
        x = rng.standard_normal(space.shape)
        live = space.flat_on(x, set(labels))
        slot = P.slot_order_flat(space, x, labels, stated)
        rev = P.slot_order_flat(space, x, labels, stated, reverse=True)
        equal_slot += int(live.tobytes() == slot.tobytes())
        differs_reversed += int(rev.tobytes() != slot.tobytes())
    return equal_slot, differs_reversed


def test_flat_in_slot_order() -> None:
    _flat_in_slot_order_positive()
    _flat_in_slot_order_signed()


def _flat_in_slot_order_positive() -> None:
    V = P.triple_space()
    equal_slot, differs = _walk_rows(V, ("angular", "spatial"), {"angular": P.ANGULAR_W, "spatial": P.SPATIAL_W})
    _require(equal_slot == N_DRAWS, f"flat_on equals the slot-order spelling on {equal_slot} of {N_DRAWS} draws")
    _require(differs >= MIN_REVERSED_DIFFERS, f"the reversed spelling differs on {differs} of {N_DRAWS} draws: the instrument does not see slot order")


def _flat_in_slot_order_signed() -> None:
    V = P.signed_pair_space()
    equal_slot, differs = _walk_rows(V, ("a", "b"), {"a": P.SIGNED_A, "b": P.SIGNED_B})
    _require(equal_slot == N_DRAWS, f"flat_on equals the slot-order spelling on {equal_slot} of {N_DRAWS} draws (signed blocks)")
    _require(differs >= MIN_REVERSED_DIFFERS, f"the reversed walk differs on {differs} of {N_DRAWS} draws on the signed pair")


def test_the_formless_slot_is_skipped_by_flat_on() -> None:
    V = P.triple_space()
    x = np.random.default_rng(3).standard_normal(V.shape)
    _require(V.flat_on(x, {"angular", "energy", "spatial"}).tobytes() == V.flat_on(x, {"angular", "spatial"}).tobytes(), "a NoForm slot is not skipped by the walk")


def _energy_window_space() -> object:
    TS = P.tensor_space()
    return TS.of_axes(P.formless_factor(P.formless_axis(3, "energy")), P.density_factor(P.diag_axis(P.SPATIAL_W, "spatial")))


def test_energy_window_is_today() -> None:
    V = _energy_window_space()
    q = np.random.default_rng(7).standard_normal(V.shape)
    today = float(np.sqrt(np.sum(q * V.apply_metric(q))))
    _require(np.float64(V.norm(q)).tobytes() == np.float64(today).tobytes(), f"norm(q) = {V.norm(q)!r} is not sqrt(sum(q * apply_metric(q))) = {today!r}")
    _require(np.float64(V.euclidean().norm(q)).tobytes() == np.float64(today).tobytes(), "euclidean().norm(q) differs from today's norm (the energy is counting, the spatial form kept)")
    counting = float(np.sqrt(np.sum(q * q)))
    _require(np.float64(V.counting().norm(q)).tobytes() == np.float64(counting).tobytes(), "counting().norm(q) is not sqrt(sum(q * q))")
    _require(counting != today, "the vector does not separate counting() from the spatial form (the fixture is inert)")
