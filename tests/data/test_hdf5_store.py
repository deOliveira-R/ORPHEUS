r"""The HDF5 store round-trips every field of an :class:`Isotope` — and refuses a stale layout.

`[M]` 2026-09-03: before #426 step 1 the serialization layer every real-library recipe
reads (`orpheus/data/micro_xs/hdf5_io.py`) had **0** test sites. Step 1 changed its
schema — ``sig2`` became the per-Legendre-order group ``sig2/L{j}`` and ``sigS`` grew
from three orders to every order the tape stores — so the round-trip is pinned here.

**Claim kind.** THEOREM about the storage contract: ``load ∘ save = id`` on every
field. COO round-trip is EXACT, so every row is ``array_equal``; a tolerance here
would be a red flag. The refusal rows pin the OTHER half of the contract: the store is
untracked and every checkout regenerates it, so a pre-#426 file must fail at the loader
with the regeneration command — never be served, never turn into a skip.

Fixtures: BE009 (8195 (n,2n) entries, 6 σ₀ columns) and H_001 (no MT=16 → the zero P0
block; 1 σ₀ column). `[M]` one ``load_isotope`` = 0.02–0.07 s.
"""
from __future__ import annotations

import h5py
import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.micro_xs import load_isotope
from orpheus.data.micro_xs.hdf5_io import H5_FORMAT, load_isotope_h5, save_isotope

pytestmark = pytest.mark.foundation

_ISOTOPES = ("BE009", "H_001")


def _write_and_read(iso, tmp_path):
    path = tmp_path / f"{iso.name}.h5"
    with h5py.File(path, "w") as f:
        save_isotope(iso, f)
    return path, load_isotope_h5(path, int(round(iso.temp)))


def _sparse_equal(a: csr_matrix, b: csr_matrix) -> bool:
    """Exact equality, pattern-insensitive (an explicit zero must not fail it)."""
    return np.array_equal(np.asarray(a.todense()), np.asarray(b.todense()))


@pytest.fixture(scope="module", params=_ISOTOPES)
def pair(request, tmp_path_factory):
    iso = load_isotope(request.param, 294)
    path, back = _write_and_read(iso, tmp_path_factory.mktemp("h5rt"))
    return iso, back, path


class TestRoundTrip:
    def test_scalar_and_vector_fields_survive(self, pair):
        iso, back, _ = pair
        assert back.aw == iso.aw
        assert back.temp == iso.temp
        for name in ("eg", "sig0", "sigC", "sigL", "sigF", "sigT", "nubar", "chi"):
            np.testing.assert_array_equal(
                np.asarray(getattr(back, name)), np.asarray(getattr(iso, name)),
                err_msg=f"{name} did not survive the HDF5 round trip",
            )

    def test_the_n2n_stack_survives_order_by_order(self, pair):
        iso, back, _ = pair
        assert len(back.sig2) == len(iso.sig2), (
            f"the loader derives the (n,2n) order count from the sig2/L{{j}} keys; "
            f"it read {len(back.sig2)}, stored {len(iso.sig2)}"
        )
        for l, (got, want) in enumerate(zip(back.sig2, iso.sig2)):
            assert _sparse_equal(got, want), f"sig2[{l}] did not round-trip"

    def test_the_scattering_stack_survives_per_order_per_sigma0(self, pair):
        iso, back, _ = pair
        assert len(back.sigS) == len(iso.sigS), (
            f"the loader derives the Legendre count from the h5 keys; it read "
            f"{len(back.sigS)}, stored {len(iso.sigS)}"
        )
        for j, (got_order, want_order) in enumerate(zip(back.sigS, iso.sigS)):
            assert len(got_order) == len(want_order), (
                f"sigS[{j}] σ₀-column count moved: {len(got_order)} vs {len(want_order)}"
            )
            for k, (got, want) in enumerate(zip(got_order, want_order)):
                assert _sparse_equal(got, want), f"sigS[{j}][{k}] did not round-trip"

    def test_the_order_axis_is_not_transposed(self, pair):
        r"""``L{j}_S{k}`` packs TWO indices; a transposed write is invisible only when
        the two counts coincide. Asserted, not skipped: `[M]` on every shipped isotope
        they differ (7 orders × {1, 6, 10} σ₀), so the row is a real discriminator, and
        the day a fixture makes it vacuous the suite must say which claim it lost."""
        iso, back, _ = pair
        n_l, n_s = len(iso.sigS), len(iso.sigS[0])
        if n_l == n_s:
            pytest.fail(
                f"fixture went vacuous: n_legendre == n_sig0 == {n_l}, so a transposed "
                f"L{{j}}_S{{k}} key is unspellable — pick an isotope whose counts differ"
            )
        assert (len(back.sigS), len(back.sigS[0])) == (n_l, n_s)

    def test_the_store_declares_its_format(self, pair):
        _, _, path = pair
        with h5py.File(path, "r") as f:
            assert f.attrs["orpheus_h5_format"] == H5_FORMAT


class TestAStaleStoreIsRefused:
    """A pre-#426 file (format 1: one ``sig2`` triplet, three orders) must fail LOUDLY.

    `[M]` the store is ``.gitignore``d and regenerated per checkout (5–12 min); recipe
    tests used to SKIP when it was missing. A schema change served silently would hand
    every recipe the wrong shape; a skip would hide the checkout that never regenerated.
    """

    def test_a_file_without_the_format_attribute_is_format_one(self, tmp_path):
        path = tmp_path / "stale.h5"
        with h5py.File(path, "w") as f:
            f.require_group("294K")
        with pytest.raises(ValueError, match=r"format 1.*convert_gxs_to_hdf5\.py"):
            load_isotope_h5(path, 294)

    def test_a_foreign_format_is_refused_with_the_regeneration_command(self, tmp_path):
        path = tmp_path / "future.h5"
        with h5py.File(path, "w") as f:
            f.attrs["orpheus_h5_format"] = H5_FORMAT + 1
            f.require_group("294K")
        with pytest.raises(ValueError, match=rf"format {H5_FORMAT + 1}.*needs format {H5_FORMAT}"):
            load_isotope_h5(path, 294)
