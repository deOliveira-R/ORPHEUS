"""HDF5 serialization for Isotope cross section data.

Each isotope+temperature is stored as an HDF5 group inside a per-element
file.  For example, ``H_001.h5`` contains groups ``/294K``, ``/350K``, etc.

The file carries ``attrs["orpheus_h5_format"] = 2``. Format 1 (until #426
step 1) stored ONE ``sig2`` triplet and three ``sigS`` orders; format 2 stores
every Legendre order the tape carries and the (n,2n) channel as a per-order
stack. A format-1 file is REFUSED by :func:`load_isotope_h5` with the
regeneration command — never silently read, never skipped.

Layout inside each temperature group::

    /{temp_K}K/
        aw          : scalar
        temp        : scalar
        eg          : (NG+1,)
        sig0        : (n_sig0,)
        sigC        : (n_sig0, NG)
        sigL        : (n_sig0, NG)
        sigF        : (n_sig0, NG)
        sigT        : (n_sig0, NG)
        nubar       : (NG,)
        chi         : (NG,)
        sig2/
            L{j}/
                row  : (nnz,)
                col  : (nnz,)
                data : (nnz,)
        sigS/
            L{j}_S{k}/
                row  : (nnz,)
                col  : (nnz,)
                data : (nnz,)
"""

from __future__ import annotations

from pathlib import Path
from typing import cast

import h5py
import numpy as np
from scipy.sparse import csr_matrix, coo_matrix

from .isotope import NG, Isotope

#: The store's format. Bumped when the layout changes so a stale cache fails
#: at the loader with the regeneration command, instead of serving a shape the
#: reader no longer expects (the store is untracked; every checkout rebuilds it).
H5_FORMAT = 2
_REGENERATE_HINT = (
    "regenerate the store with "
    "`.venv/bin/python orpheus/data/micro_xs/convert_gxs_to_hdf5.py`"
)


def _order_key(order: int, sigma0: int | None = None) -> str:
    """The HDF5 group name of one Legendre block: ``L{j}`` or ``L{j}_S{k}``."""
    return f"L{order}" if sigma0 is None else f"L{order}_S{sigma0}"


def _n_orders(group: h5py.Group) -> int:
    """How many Legendre orders a block group stores, read off its ``L{j}…`` keys."""
    return max(int(k.split("_")[0][1:]) for k in group.keys()) + 1


def save_isotope(iso: Isotope, h5file: h5py.File) -> None:
    """Write one Isotope into an open HDF5 file as a temperature group."""
    h5file.attrs["orpheus_h5_format"] = H5_FORMAT
    temp_K = int(round(iso.temp))
    grp = h5file.require_group(f"{temp_K}K")

    grp.attrs["aw"] = iso.aw
    grp.attrs["temp"] = iso.temp

    for name, data in [
        ("eg", iso.eg), ("sig0", iso.sig0),
        ("sigC", iso.sigC), ("sigL", iso.sigL),
        ("sigF", iso.sigF), ("sigT", iso.sigT),
        ("nubar", iso.nubar), ("chi", iso.chi),
    ]:
        if name in grp:
            del grp[name]
        grp.create_dataset(name, data=data, compression="gzip", compression_opts=4)

    # Sparse matrices: store as COO triplets, one group per Legendre order
    sig2_grp = grp.require_group("sig2")
    for j, mat in enumerate(iso.sig2):
        _save_sparse(sig2_grp, _order_key(j), mat)

    sig_s_grp = grp.require_group("sigS")
    for j, legendre_mats in enumerate(iso.sigS):
        for k, mat in enumerate(legendre_mats):
            _save_sparse(sig_s_grp, _order_key(j, k), mat)


def load_isotope_h5(path: Path, temp_K: int) -> Isotope:
    """Load an Isotope from an HDF5 file for a given temperature."""
    with h5py.File(path, "r") as f:
        stored = int(f.attrs.get("orpheus_h5_format", 1))
        if stored != H5_FORMAT:
            raise ValueError(
                f"{path.name} is HDF5 store format {stored}; this reader needs "
                f"format {H5_FORMAT} (per-Legendre-order (n,2n) stack, every "
                f"tape order kept — #426) — {_REGENERATE_HINT}"
            )
        grp = cast(h5py.Group, f[f"{temp_K}K"])

        aw = float(cast(np.floating, grp.attrs["aw"]))
        temp = float(cast(np.floating, grp.attrs["temp"]))
        eg = cast(h5py.Dataset, grp["eg"])[:]
        ng = len(eg) - 1
        sig0 = cast(h5py.Dataset, grp["sig0"])[:]
        n_sig0 = len(sig0)

        sigC = cast(h5py.Dataset, grp["sigC"])[:]
        sigL = cast(h5py.Dataset, grp["sigL"])[:]
        sigF = cast(h5py.Dataset, grp["sigF"])[:]
        sigT = cast(h5py.Dataset, grp["sigT"])[:]
        nubar = cast(h5py.Dataset, grp["nubar"])[:]
        chi = cast(h5py.Dataset, grp["chi"])[:]

        # Reconstruct sig2: [legendre] — the order count is read off the keys,
        # exactly as sigS's is below, so the store never needs to state it.
        sig2_grp = cast(h5py.Group, grp["sig2"])
        sig2 = [_load_sparse(sig2_grp, _order_key(j), ng=ng) for j in range(_n_orders(sig2_grp))]

        # Reconstruct sigS structure: [legendre][sig0_idx]
        sig_s_grp = cast(h5py.Group, grp["sigS"])
        sigS = [
            [_load_sparse(sig_s_grp, _order_key(j, k), ng=ng) for k in range(n_sig0)]
            for j in range(_n_orders(sig_s_grp))
        ]

    name = path.stem + f"_{temp_K}K"
    return Isotope(
        name=name, aw=aw, temp=temp, eg=eg, sig0=sig0,
        sigC=sigC, sigL=sigL, sigF=sigF, sigT=sigT,
        nubar=nubar, chi=chi, sigS=sigS, sig2=sig2,
    )


def _save_sparse(parent: h5py.Group, name: str, mat: csr_matrix) -> None:
    """Save a sparse matrix as COO triplets."""
    grp = parent.require_group(name)
    coo = mat.tocoo()
    for key in ("row", "col", "data"):
        if key in grp:
            del grp[key]
    grp.create_dataset("row", data=coo.row.astype(np.int32), compression="gzip")
    grp.create_dataset("col", data=coo.col.astype(np.int32), compression="gzip")
    grp.create_dataset("data", data=coo.data, compression="gzip")


def _load_sparse(parent: h5py.Group, name: str, ng: int = NG) -> csr_matrix:
    """Load a sparse matrix from COO triplets."""
    grp = cast(h5py.Group, parent[name])
    row = cast(h5py.Dataset, grp["row"])[:]
    col = cast(h5py.Dataset, grp["col"])[:]
    data = cast(h5py.Dataset, grp["data"])[:]
    # coo_matrix(...).tocsr() returns a csr_matrix at runtime (matrix lineage);
    # pyright infers the untyped tocsr() body as csr_array, so cast to the truth.
    return cast(csr_matrix, coo_matrix((data, (row, col)), shape=(ng, ng)).tocsr())
