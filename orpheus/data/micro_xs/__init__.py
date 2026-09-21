"""Microscopic cross section data loading.

Loads isotope data from the HDF5 store, a derived cache of the tracked GENDF
(``.GXS``) tapes that ``convert_gxs_to_hdf5.py`` writes and that the loader
materialises file by file on first use.
"""

from pathlib import Path

from .isotope import Isotope
from .hdf5_io import load_isotope_h5
from .convert_gxs_to_hdf5 import convert_one

_HDF5_DIR = Path(__file__).resolve().parent


def load_isotope(name: str, temp_K: int) -> Isotope:
    """Load an isotope from the HDF5 store, materialising its file from the
    tracked ``{name}.GXS`` tape when the store does not carry it yet.

    The store is untracked (``.gitignore``: ``*.h5``), a cache of the tapes,
    so a fresh checkout — a CI runner, a new machine — builds each file the
    first time an isotope is asked for (`[M]` 2026-09-21: about 1.1 s per MB
    of tape, O_016 in 32 s) and every later load reads it in milliseconds. A
    file of a stale format is not rebuilt here: it refuses at the reader with
    the regeneration command, the contract ``tests/data/test_hdf5_store.py``
    pins.
    """
    h5_path = _HDF5_DIR / f"{name}.h5"
    if not h5_path.exists():
        convert_one(name, _HDF5_DIR)
    return load_isotope_h5(h5_path, temp_K)


__all__ = ["Isotope", "load_isotope"]
