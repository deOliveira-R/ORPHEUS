r"""The lossless ingest (#426 step 1) APPENDED Legendre orders and moved no stored value.

A digest ledger of the shipped store as it stood before step 1 — captured 2026-09-03 on
``main`` ``1e02f6b1`` by ``scratch/_426_capture_ingest_ledger.py`` from the format-1 h5
files — pins, for all 13 isotopes at 294 K: ``sig0``, ``sigT``, ``chi``, ``nubar``, the
row sums and nnz of ``sigS[l][k]`` for the three orders that existed, and the row sum and
nnz of the P0 (n,2n) block. Every digest must reproduce from the regenerated store, and
every stack must be at least as long as it was.

**Claim kind.** A frozen-reference regression pin (`vv` ladder: the reference is the tree's
own earlier reading, so this proves *unchanged*, not *correct*; correctness of the values
is the tape-physics and yield rows in ``test_n2n_yield_convention.py``). It exists because
the h5 store is untracked: a checkout that regenerates it must land on the same P0..P2
values the whole verification history was measured against.
"""
from __future__ import annotations

import hashlib
import json

import numpy as np
import pytest

from orpheus.data.micro_xs import load_isotope

from . import DATA_TREE_DATA

pytestmark = pytest.mark.foundation

_LEDGER = DATA_TREE_DATA / "pre_426_ingest_ledger.json"


def _digest(a) -> dict:
    a = np.ascontiguousarray(np.asarray(a, dtype=float))
    return {"sha256": hashlib.sha256(a.tobytes()).hexdigest(), "shape": list(a.shape)}


@pytest.fixture(scope="module")
def ledger() -> dict:
    with open(_LEDGER) as f:
        return json.load(f)


def _isotopes(ledger: dict) -> list[str]:
    return sorted(k for k in ledger if not k.startswith("_"))


@pytest.fixture(scope="module")
def isotopes(ledger):
    return {name: load_isotope(name, 294) for name in _isotopes(ledger)}


@pytest.mark.parametrize("name", [
    "BE009", "B_010", "B_011", "H_001", "NA023", "O_016", "U_235", "U_238",
    "ZR090", "ZR091", "ZR092", "ZR094", "ZR096",
])
class TestTheStoredValuesDidNotMove:
    def test_vectors(self, ledger, isotopes, name):
        iso, want = isotopes[name], ledger[name]
        for field in ("sig0", "sigT", "chi", "nubar"):
            assert _digest(getattr(iso, field)) == want[field], f"{name}.{field} moved"

    def test_scattering_orders_zero_to_two_at_every_sigma0(self, ledger, isotopes, name):
        iso, want = isotopes[name], ledger[name]
        assert len(iso.sigS) >= want["n_legendre"], f"{name}: the stack got SHORTER"
        assert iso.n_sig0 == want["n_sig0"]
        for key, digest in want["sigS_rowsum"].items():
            l, k = (int(t[1:]) for t in key.split("_"))
            assert _digest(np.asarray(iso.sigS[l][k].sum(axis=1)).ravel()) == digest, (
                f"{name}.sigS[{l}][{k}] row sums moved"
            )
            assert iso.sigS[l][k].nnz == want["sigS_nnz"][key], f"{name}.sigS[{l}][{k}] nnz moved"

    def test_the_p0_n2n_block(self, ledger, isotopes, name):
        iso, want = isotopes[name], ledger[name]
        assert _digest(np.asarray(iso.sig2[0].sum(axis=1)).ravel()) == want["sig2_rowsum"], (
            f"{name}.sig2[0] row sums moved"
        )
        assert iso.sig2[0].nnz == want["sig2_nnz"]


class TestTheOrdersWereAppended:
    """What step 1 ADDED: every scattering section's NL, and the (n,2n) stack, are kept."""

    def test_every_isotope_stores_seven_scattering_orders(self, isotopes):
        """`[M]` NL = 7 on every shipped MT=2 section (13 of 13 tapes)."""
        short = {n: len(i.sigS) for n, i in isotopes.items() if len(i.sigS) != 7}
        assert not short, f"expected 7 stored scattering orders everywhere; got {short}"

    def test_the_n2n_stack_has_the_tapes_order(self, isotopes):
        """`[M]` MT=16: NL = 7 on 10 tapes, NL = 1 on NA023, absent on B_010 and H_001
        (the zero P0 block) — the ingest invents nothing and pads nothing."""
        got = {n: len(i.sig2) for n, i in isotopes.items()}
        want = {n: 7 for n in got}
        want.update({"NA023": 1, "B_010": 1, "H_001": 1})
        assert got == want

    def test_a_channel_the_tape_lacks_is_the_zero_p0_block(self, isotopes):
        for name in ("B_010", "H_001"):
            assert isotopes[name].sig2[0].nnz == 0

    def test_a_mixture_of_unequal_stacks_pads_the_short_one(self, isotopes):
        """Σ_i N_i σ_{i,ℓ} order by order; an isotope with no ℓ ≥ 1 (n,2n) contributes zero there.

        Be-9 (7 orders) + B-10 (the zero P0 block, 1 order) at equal densities: the
        mixture's stack has 7 orders and every block is exactly Be-9's — B-10 adds
        nothing at any ℓ, and the sum does not fail on the short stack (ruling O-1/O-2:
        pad, never clamp). `[M]` 2026-09-03: a `_macroscopic_stack` that skipped the
        padding would raise or drop orders here.
        """
        import contextlib, io
        from orpheus.data.macro_xs.mixture import compute_macro_xs

        be, b10 = isotopes["BE009"], isotopes["B_010"]
        with contextlib.redirect_stdout(io.StringIO()):
            mixed = compute_macro_xs([be, b10], np.array([1.0, 1.0]))
            alone = compute_macro_xs([be], np.array([1.0]))
        assert len(mixed.Sig2) == 7 == len(alone.Sig2)
        for l, (got, want) in enumerate(zip(mixed.Sig2, alone.Sig2)):
            assert np.array_equal(got.todense(), want.todense()), f"Sig2[{l}] moved when B-10 joined"
        assert len(mixed.SigS) == 7, "the scattering stack keeps every order too"
