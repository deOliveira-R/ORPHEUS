r"""Extractors: :class:`Mixture` → numpy arrays for legacy consumers.

Production solvers (CP, SN, MOC) consume :class:`Mixture` directly.
Semi-analytical reference solvers (F_N, PS-1982 wrapper, transfer-
matrix k_inf) were written before the protocol and accept raw numpy
arrays of :math:`\\Sigma_t`, :math:`\\Sigma_s`,
:math:`\\nu \\Sigma_f`, :math:`\\chi`. This module bridges the gap.
"""
from __future__ import annotations

import numpy as np

from orpheus.data.macro_xs.mixture import Mixture


def mixture_to_fn_arrays(
    mixture: Mixture,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    r"""Extract ``(sigma_t, sigma_s_p0, nu_sigma_f, chi)`` from a Mixture.

    The semi-analytical reference solvers (F_N method, PS-1982 wrapper,
    :func:`orpheus.derivations.common.eigenvalue.kinf_homogeneous`)
    consume:

    * ``sigma_t``: dense ``(NG,)`` total XS.
    * ``sigma_s_p0``: dense ``(NG, NG)`` P_0 scattering matrix in
      ORPHEUS ``[from, to]`` convention.
    * ``nu_sigma_f``: dense ``(NG,)`` production XS.
    * ``chi``: dense ``(NG,)`` fission spectrum.

    Parameters
    ----------
    mixture : Mixture
        Production-protocol Mixture (typically from
        :func:`orpheus.derivations.common.xs_library.make_mixture`).

    Returns
    -------
    sigma_t, sigma_s_p0, nu_sigma_f, chi : np.ndarray
        Numpy arrays with the shapes documented above. Returns dense
        arrays — callers may ``.copy()`` if they need to mutate.

    Notes
    -----
    ``Mixture.SigS`` is a list of sparse matrices indexed by Legendre
    order; this extractor returns ``SigS[0].toarray()`` (the P_0
    moment). Higher Legendre moments are out of scope for the
    isotropic-scattering first-slice cases.
    """
    sigma_t = np.asarray(mixture.SigT, dtype=float)
    if not mixture.SigS:
        raise ValueError("Mixture has no scattering matrices (SigS is empty)")
    sigma_s_p0 = mixture.SigS[0].toarray().astype(float)
    nu_sigma_f = np.asarray(mixture.SigP, dtype=float)
    chi = np.asarray(mixture.chi, dtype=float)
    return sigma_t, sigma_s_p0, nu_sigma_f, chi
