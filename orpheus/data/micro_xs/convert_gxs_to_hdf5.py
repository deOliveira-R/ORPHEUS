#!/usr/bin/env python3
"""Convert all GENDF (.GXS) files to HDF5 format.

Reads from ``01.Micro.XS.421g/*.GXS`` and writes one ``.h5`` file per
element into ``data/micro_xs/``.  Each HDF5 file contains all temperatures.
"""

from pathlib import Path

import h5py

from orpheus.data.micro_xs.gendf import convert_gxs, _GXS_DIR
from orpheus.data.micro_xs.hdf5_io import save_isotope

OUTPUT_DIR = Path(__file__).resolve().parent  # same directory as this script


def convert_one(name: str, out_dir: Path = OUTPUT_DIR) -> Path:
    """Materialise ``{name}.h5`` under ``out_dir`` from the tracked ``{name}.GXS``
    tape: every temperature the tape carries, in the store's current format.

    The ONE producer of a store file — the whole-library rebuild below and the
    loader's first-use materialisation (:func:`orpheus.data.micro_xs.load_isotope`)
    both call it, so the two paths cannot drift. `[M]` 2026-09-21: about 1.1 s
    per MB of tape on the maintainer's machine (O_016 32 s, U_235 74 s).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    h5_path = out_dir / f"{name}.h5"

    print(f"Converting {name}:")
    isotopes = convert_gxs(name)

    with h5py.File(h5_path, "w") as f:
        f.attrs["element"] = name
        f.attrs["n_temperatures"] = len(isotopes)
        for iso in isotopes:
            save_isotope(iso, f)
            temp_K = int(round(iso.temp))
            print(f"    {temp_K}K written")

    size_mb = h5_path.stat().st_size / 1e6
    print(f"  -> {h5_path.name} ({size_mb:.1f} MB, {len(isotopes)} temperatures)\n")
    return h5_path


def main():
    gxs_files = sorted(_GXS_DIR.glob("*.GXS"))
    print(f"Found {len(gxs_files)} GXS files in {_GXS_DIR}\n")

    for gxs_path in gxs_files:
        convert_one(gxs_path.stem, OUTPUT_DIR)  # e.g. "H_001", "U_235"

    print("Done. All HDF5 files in:", OUTPUT_DIR.resolve())


if __name__ == "__main__":
    main()
