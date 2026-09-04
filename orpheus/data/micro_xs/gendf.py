"""GENDF (GXS) file parser.

Reads the IAEA 421-group GENDF files directly — no CSV intermediary needed.
The format is fixed-width 80-column ENDF-style records:

    Columns  1-66 : 6 data fields, 11 characters each
    Columns 67-70 : MAT number
    Columns 71-72 : MF  number
    Columns 73-75 : MT  number
    Columns 76-80 : line sequence number

Data fields use a compact float notation where 'E' is omitted,
e.g. ``1.001000+3`` means ``1.001E+3``.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import replace
from pathlib import Path
from typing import cast

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix, diags

from .isotope import NG, Isotope

_GXS_DIR = Path(__file__).resolve().parent


# ---------------------------------------------------------------------------
# Low-level GENDF parsing
# ---------------------------------------------------------------------------

def _parse_gendf_field(s: str) -> float:
    """Parse one 11-character GENDF data field into a float.

    Handles the compact notation where 'E' is omitted:
    ``' 1.001000+3'`` → ``1.001e+3``.
    """
    s = s.strip()
    if not s:
        return 0.0
    # Insert 'E' before +/- sign that follows a digit (but not at position 0)
    s = re.sub(r"(\d)([+-])", r"\1E\2", s)
    return float(s)


def _parse_gendf(path: Path) -> np.ndarray:
    """Read a GXS file into a numeric matrix (n_lines, 10).

    Columns: [data1..data6, MAT, MF, MT, line_seq].
    Mirrors what MATLAB's ``importdata('file.CSV', ';')`` produces.
    """
    rows = []
    with open(path) as f:
        for line in f:
            if len(line.rstrip("\n")) < 75:
                continue
            data = [_parse_gendf_field(line[i * 11 : (i + 1) * 11]) for i in range(6)]
            mat = int(line[66:70])
            mf = int(line[70:72])
            mt = int(line[72:75])
            seq = int(line[75:80])
            rows.append(data + [mat, mf, mt, seq])
    return np.array(rows)


# ---------------------------------------------------------------------------
# Cross section extraction (ports of MATLAB extract_mf3, extract_mf6)
# ---------------------------------------------------------------------------

def _extract_n_words(n: int, i_row: int, m: np.ndarray) -> tuple[np.ndarray, int]:
    """Read *n* consecutive data words starting at row ``i_row``.

    Words are packed 6 per row in columns 0-5.  Returns the extracted
    values and the row index of the last row read.
    """
    a = []
    row = i_row
    full_lines = n // 6
    for _ in range(full_lines):
        a.extend(m[row, :6])
        row += 1
    remainder = n - full_lines * 6
    if remainder > 0:
        a.extend(m[row, :remainder])
    else:
        row -= 1  # last full row was the final one
    return np.array(a), row


def _extract_mf3(mt: int, temp_idx: int, m: np.ndarray) -> np.ndarray | None:
    """Extract MF=3 cross sections for reaction *mt* at temperature index *temp_idx*.

    Returns (n_sig0, NG) array, or None if the reaction is absent.
    """
    n_row = m.shape[0]
    n_temp = 0
    i_row_found = 0

    for i in range(n_row):
        if m[i, 7] == 3 and m[i, 8] == mt and m[i, 9] == 1:
            n_temp += 1
            if n_temp == temp_idx:
                i_row_found = i + 1
                break

    if i_row_found == 0:
        return None

    n_sig0 = int(m[i_row_found - 1, 3])
    n_lgn = int(m[i_row_found - 1, 2])
    sig = np.zeros((n_sig0, NG))

    i = i_row_found + 1
    while i < n_row and m[i, 7] == 3 and m[i, 8] == mt:
        ig = int(m[i - 1, 5]) - 1  # 0-based group index
        n_words = n_sig0 * n_lgn * 2
        a, i_new = _extract_n_words(n_words, i, m)
        # The second half contains the XS values for each sigma-zero
        sig[:, ig] = a[n_sig0 * n_lgn : n_sig0 * n_lgn + n_sig0]
        i = i_new + 2
    return sig


def _extract_mf6(
    mt: int, temp_idx: int, m: np.ndarray
) -> tuple[np.ndarray, np.ndarray, dict[tuple[int, int], np.ndarray]] | None:
    """Extract MF=6 transfer matrix for reaction *mt* at temperature index *temp_idx*.

    Returns (ifrom, ito, sig_dict) where sig_dict[(legendre, sig0_idx)] is
    a 1-D array of non-zero values.  ifrom/ito are 1-based group indices
    (matching MATLAB convention for later sparse matrix construction).
    Returns None if the reaction is absent.

    Every Legendre order the section stores is returned — ``NL`` keys per
    sigma-zero column — and nothing is padded: a section with ``NL = 1`` is
    the evaluation declaring isotropy, and its higher moments are exactly
    zero for every consumer, which pads (#426 step 1; until then ``NL = 1``
    was padded to three zero keys here and every channel was cut to three
    orders downstream).
    """
    n_row = m.shape[0]
    i = 0
    n_temp = 0
    ifrom_list: list[int] = []
    ito_list: list[int] = []
    sig: dict[tuple[int, int], list[float]] = {}

    # n_lgn/n_sig0 are set from the section header (the SEQ==1 record), which
    # also increments n_temp.  The read sites below are all guarded by
    # ``n_temp == temp_idx`` with temp_idx >= 1, so a read can only happen
    # after the header bound both — these defaults are provably never read and
    # exist only so the control flow is statically sound (no UnboundLocalError
    # on a malformed file, and no false "possibly unbound" from the checker).
    n_lgn = 0
    n_sig0 = 0

    while i < n_row and m[i, 6] != -1:
        if m[i, 7] == 6 and m[i, 8] == mt:
            if m[i, 9] == 1:  # first record of this MF/MT section
                n_lgn = int(m[i, 2])
                n_sig0 = int(m[i, 3])
                i += 1
                n_temp += 1

            ng2 = int(m[i, 2])      # number of secondary positions
            ig2lo = int(m[i, 3])    # lowest nonzero group (1-based)
            nw = int(m[i, 4])       # words to read
            ig = int(m[i, 5])       # current group (1-based)

            i += 1
            a, i_new = _extract_n_words(nw, i, m)
            i = i_new

            if n_temp == temp_idx:
                k = n_lgn * n_sig0  # skip flux words
                for i_to in range(ig2lo, ig2lo + ng2 - 1):
                    ifrom_list.append(ig)
                    ito_list.append(i_to)
                    for i_sig0 in range(n_sig0):
                        for i_lgn in range(n_lgn):
                            k += 1
                            sig.setdefault((i_lgn, i_sig0), []).append(a[k - 1])
        i += 1

    if n_temp == 0:
        return None

    ifrom = np.array(ifrom_list)
    ito = np.array(ito_list)
    sig_arrays = {key: np.array(vals) for key, vals in sig.items()}
    return ifrom, ito, sig_arrays


def _legendre_order(data: dict[tuple[int, int], np.ndarray]) -> int:
    """The number of Legendre orders a MF=6 section stores (its ``NL``)."""
    return 1 + max(legendre for legendre, _ in data)


def _strip_transfer_yield(
    transfer_stack: Sequence[csr_matrix], reaction_xs: np.ndarray | None, *, mt: int
) -> list[csr_matrix]:
    r"""Renormalise a MF=6 transfer STACK onto its MF=3 reaction XS.

    GENDF's MF=6 records hold :math:`\sigma(E)\,y(E)\,f(E\to E')` — the
    transfer matrix carries the reaction's **yield** — while MF=3 holds the
    un-multiplied group cross section :math:`\sigma`. For a yield-1 channel
    (elastic, inelastic) the two agree and this is the identity; for
    :math:`(n,2n)` the MF=6 row sum is :math:`2\sigma` and the multiplicity
    must come off before the matrix reaches a consumer that applies it
    itself.

    Scaling each row onto ``reaction_xs`` divides the yield out **without
    naming its value**, which is the point: the multiplicity is a physics
    constant owned once, by
    :attr:`~orpheus.transport.kernels.N2NKernel.multiplicity`, in a package
    this one may not import. It also makes the channel's reaction rate
    exactly consistent with the MF=3 tabulation every other channel's cross
    section is read from, rather than merely consistent to the file's
    six-digit rounding.

    The stack carries ONE yield. GENDF stores :math:`\sigma_\ell =
    \sigma(E)\,y(E)\,f_\ell(E\to E')` for every :math:`\ell` (ENDF-102
    Eq. (6.1)/(6.3) with NJOY Eq. (242)): the per-row scale is derived from
    :math:`\ell = 0` — whose row sum is :math:`y\,\sigma` — and applied to
    every order. An :math:`\ell \ge 1` row sums to :math:`y\,\sigma\,
    \langle P_\ell \rangle` and could not be normalised on its own; the
    integer-yield admission below is therefore a statement about the P0
    block alone.

    Raises
    ------
    ValueError
        If ``reaction_xs`` is absent while a transfer matrix exists (a
        MF=6 section with no MF=3 partner is a malformed tape — `[M]` all
        11 shipped files carrying MT=16 have both), if the reaction XS
        vanishes where the transfer matrix does not, or if the aggregate
        yield is not a positive integer. The last is the admission
        contract: a yield of 1.5 means this reader has misunderstood the
        tape's convention, and silently renormalising would bury that.
        A real ``raise`` rather than an ``assert`` — the canonical runner
        is ``python -O``, which strips ``assert`` at compile time.
    """
    p0 = transfer_stack[0]
    rows = np.asarray(p0.sum(axis=1)).ravel()
    live = rows > 0.0
    if not live.any():
        return list(transfer_stack)
    if reaction_xs is None:
        raise ValueError(
            f"GENDF MT={mt}: a MF=6 transfer matrix is present but its "
            f"MF=3 cross section is absent, so the record's yield cannot "
            f"be divided out."
        )
    sigma = np.asarray(reaction_xs)[0]
    if np.any(sigma[live] <= 0.0):
        bad = int(np.flatnonzero(live & (sigma <= 0.0))[0])
        raise ValueError(
            f"GENDF MT={mt}: the MF=6 transfer matrix is non-zero in group "
            f"{bad} where the MF=3 cross section is {sigma[bad]:g}."
        )
    aggregate = float(rows[live].sum() / sigma[live].sum())
    nearest = round(aggregate)
    if nearest < 1 or abs(aggregate - nearest) > 1.0e-3:
        raise ValueError(
            f"GENDF MT={mt}: rowsum(MF=6)/sigma(MF=3) aggregates to "
            f"{aggregate:.6f}, which is not a positive integer yield. "
            f"This reader assumes MF=6 carries sigma*y*f; a non-integral "
            f"ratio means that assumption does not hold for this tape."
        )
    scale = np.ones(NG)
    scale[live] = sigma[live] / rows[live]
    yield_inverse = diags(scale)
    return [cast(csr_matrix, yield_inverse @ moment) for moment in transfer_stack]


# ---------------------------------------------------------------------------
# High-level: GXS → Isotope
# ---------------------------------------------------------------------------

_IG_THRESH = 95  # last group of thermal energy (E ≈ 4 eV)

#: Added to every stored MF=6 scattering value so an exactly-zero entry survives
#: coo→csr canonicalisation. This is what makes all ``n_sig0`` columns of ONE
#: Legendre order share one sparsity pattern — the assumption ``interp_sig_s``
#: reads positions under (it takes ``sigS[l][0]``'s coordinates and interpolates
#: every column at those positions). Load-bearing at every ℓ, not only P0: the
#: higher orders have genuine exact zeros. ⚠ It does NOT make the orders agree
#: with each other — `[M]` NA023's MT=91 section stores NL = 1, so its L0 carries
#: 16809 nnz against 13196 at ℓ ≥ 1; the other 12 isotopes agree only because
#: every scattering section on their tapes stores NL = 7. (MATLAB's value.)
_SPARSITY_EPSILON = 1e-30


def convert_gxs(name: str) -> list[Isotope]:
    """Convert a GXS file to a list of Isotope objects (one per temperature).

    Parameters
    ----------
    name : str
        Isotope identifier matching the GXS filename, e.g. ``"H_001"``,
        ``"U_235"``, ``"ZR090"``.

    Returns
    -------
    list[Isotope] — one entry per temperature found in the file.
    """
    path = _GXS_DIR / f"{name}.GXS"
    if not path.exists():
        raise FileNotFoundError(f"No GXS file: {path}")

    print(f"  Parsing {path.name}...", end=" ", flush=True)
    m = _parse_gendf(path)
    print(f"{m.shape[0]} records.", flush=True)

    # --- Header: temperatures, sigma-zeros, energy grid (MF=1, MT=451) ---
    temps: list[float] = []
    for i in range(m.shape[0]):
        if m[i, 7] == 1 and m[i, 8] == 451 and m[i, 9] == 2:
            temps.append(m[i, 0])

    n_sig0 = int(m[1, 3])
    header_words = 1 + n_sig0 + (NG + 1)
    a, _ = _extract_n_words(header_words, 3, m)
    sig0 = a[1 : 1 + n_sig0]
    eg = a[1 + n_sig0 : 1 + n_sig0 + NG + 1]
    aw = m[1, 1] * 1.008664916  # convert to amu

    isotopes = []
    for i_temp, temp in enumerate(temps, start=1):
        iso = _build_isotope(name, temp, i_temp, m, aw, eg, sig0, n_sig0)
        # Normalise NJOY thermal-first → ORPHEUS canonical fast-first at the
        # ingest boundary (see _to_canonical_group_order).
        isotopes.append(_to_canonical_group_order(iso))

    return isotopes


def _build_isotope(
    name: str,
    temp: float,
    i_temp: int,
    m: np.ndarray,
    aw: float,
    eg: np.ndarray,
    sig0: np.ndarray,
    n_sig0: int,
) -> Isotope:
    """Build one Isotope from parsed GENDF data for a given temperature."""
    # --- MF=3 reactions ---
    sigC = _extract_mf3(102, i_temp, m)  # radiative capture
    if sigC is None:
        sigC = np.zeros((n_sig0, NG))
    elif sigC.shape[0] == 1 and n_sig0 > 1:
        sigC = np.tile(sigC, (n_sig0, 1))

    sigL_raw = _extract_mf3(107, i_temp, m)  # (n,alpha)
    if sigL_raw is None:
        sigL = np.zeros((n_sig0, NG))
    elif sigL_raw.shape[0] == 1 and n_sig0 > 1:
        sigL = np.tile(sigL_raw, (n_sig0, 1))
    else:
        sigL = sigL_raw

    sigF_raw = _extract_mf3(18, i_temp, m)  # fission
    if sigF_raw is None:
        sigF = np.zeros((n_sig0, NG))
    else:
        sigF = sigF_raw

    nubar_raw = _extract_mf3(452, i_temp, m)  # nubar
    nubar = nubar_raw[0] if nubar_raw is not None else np.zeros(NG)

    # --- (n,2n) Legendre stack: MF=6, MT=16 ---
    # GENDF's MF=6 stores sigma(E)*y(E)*f(E->E'), i.e. the transfer matrix
    # carries the reaction's YIELD; MF=3 stores the un-multiplied reaction
    # cross section. For MT=16 that yield is 2, and `[M]` on all 11 shipped
    # files carrying the section the ratio rowsum(MF=6)/sigma(MF=3) is
    # 2.000000 (worst per-group departure 2.8e-2, in a threshold group where
    # sigma is ~1e-4 and the file's 6-digit fields round hardest).
    #
    # Every consumer downstream reads `sig2[l]` as the REACTION stack with no
    # multiplicity folded in — `SigT`/`absorption_xs` add the P0 row sum ONCE
    # (one neutron absorbed per event) and `N2NKernel.emission_matrix`
    # applies the factor itself (two neutrons emitted). So the yield is
    # divided out HERE, at the definition site, rather than at each of them.
    #
    # The division is spelled as a RENORMALISATION ONTO MF=3, not as a
    # literal `/ 2`: the multiplicity is a physics constant with exactly one
    # home in this tree (`N2NKernel.multiplicity`), which the data layer must
    # not import (it sits a layer up) and must not duplicate. Normalising to
    # the tabulated cross section removes whatever yield the file carries
    # without this module ever naming its value — and it makes the (n,2n)
    # reaction rate exactly consistent with the MF=3 tabulation that every
    # other channel's XS is read from. (Issue #427; the pre-fix tree fed
    # `2*Sigma_16` to a consumer set expecting `Sigma_16`, so removal was
    # doubled and the emission quadrupled.)
    n2n = _extract_mf6(16, i_temp, m)
    sig2: list[csr_matrix]
    if n2n is not None:
        ifrom2, ito2, sig2_data = n2n
        # The whole Legendre stack the section stores (NL = 7 on 10 of the 11
        # shipped files carrying MT=16, NL = 1 on NA023), one sigma-zero column:
        # a threshold channel is not self-shielded. Until #426 step 1 only
        # `(0, 0)` was kept and the emission's anisotropy was discarded here.
        # coo_matrix(...).tocsr() is a csr_matrix at runtime (matrix lineage);
        # the untyped tocsr() body makes pyright infer csr_array, so cast.
        sig2 = [
            cast(
                csr_matrix,
                coo_matrix(
                    (sig2_data[(i_lgn, 0)], (ifrom2 - 1, ito2 - 1)), shape=(NG, NG)
                ).tocsr(),
            )
            for i_lgn in range(_legendre_order(sig2_data))
        ]
        sig2 = _strip_transfer_yield(sig2, _extract_mf3(16, i_temp, m), mt=16)
    else:
        sig2 = [csr_matrix((NG, NG))]

    # --- Scattering matrices: elastic + inelastic + thermal ---
    # Every Legendre order any scattering section stores is kept (#426 step 1;
    # a hard-coded three cut elastic P3..P6 here until then). A section with
    # fewer orders than the widest contributes exactly zero above its own NL.
    elastic = _extract_mf6(2, i_temp, m)                       # MF=6, MT=2
    inelastic = [
        section
        for mt in range(51, 92)                                # MF=6, MT=51..91
        if (section := _extract_mf6(mt, i_temp, m)) is not None
    ]
    thermal_mt = 222 if name.startswith("H_001") else 221      # H-in-water / free gas
    thermal = _extract_mf6(thermal_mt, i_temp, m)
    sections = [s for s in (elastic, *inelastic, thermal) if s is not None]
    n_legendre = max((_legendre_order(s[2]) for s in sections), default=1)

    sigS = _init_scattering(elastic, n_sig0, n_legendre)
    for section in inelastic:
        _accumulate_scattering(sigS, section, n_sig0, sigma_zero_independent=True)
    if thermal is not None:
        _accumulate_scattering(sigS, thermal, n_sig0, sigma_zero_independent=True)

    # --- Fission spectrum (chi): MF=6, MT=18 ---
    chi = _extract_chi(i_temp, m)

    # --- Total cross section (computed from components) ---
    # The P0 row sum IS the (n,2n) reaction XS (one neutron absorbed per event);
    # the ℓ ≥ 1 moments integrate to zero over angle and never enter a removal.
    n2n_reaction_xs = np.array(sig2[0].sum(axis=1)).ravel()
    sigT = np.zeros((n_sig0, NG))
    for i_sig0 in range(n_sig0):
        scattering_xs = np.array(sigS[0][i_sig0].sum(axis=1)).ravel()
        sigT[i_sig0] = sigC[i_sig0] + sigF[i_sig0] + sigL[i_sig0] + scattering_xs + n2n_reaction_xs

    temp_K = int(round(temp))
    return Isotope(
        name=f"{name}_{temp_K}K",
        aw=aw,
        temp=temp,
        eg=eg,
        sig0=sig0,
        sigC=sigC,
        sigL=sigL,
        sigF=sigF,
        sigT=sigT,
        nubar=nubar,
        chi=chi,
        sigS=sigS,
        sig2=sig2,
    )


# NJOY/GENDF stores group 0 = THERMAL (energies ASCENDING); the ORPHEUS canonical
# convention is group 0 = FASTEST (highest E), ``eg`` DESCENDING, and
# ``SigS[g_from, g_to]`` downscatter in the UPPER triangle. We normalise the
# foreign NJOY order ONCE here, at the data-ingest boundary, so every downstream
# consumer (Mixture, cell_xs, every solver) is order-transparent.
# See docs/theory/foundations/cross_section_data.rst (canonical group convention).


def _reverse_groups_2d(m: csr_matrix) -> csr_matrix:
    """Reverse BOTH group axes of a square ``[g_from, g_to]`` matrix (sparse-preserving).

    Reversing both axes keeps the matrix self-consistent (row = from-group,
    col = to-group) while moving downscatter from the lower to the upper
    triangle — the canonical fast-first structure.
    """
    n = cast("tuple[int, ...]", m.shape)[0]
    perm = np.arange(n - 1, -1, -1)
    return cast(csr_matrix, m[perm][:, perm])


def _to_canonical_group_order(iso: Isotope) -> Isotope:
    """Reverse the energy-group axis: NJOY thermal-first → ORPHEUS fast-first.

    The single foreign-input normalisation at the GENDF ingest boundary. Reverses
    every group-indexed array — ``eg`` (→ descending), the ``(n_sig0, NG)`` vector
    cross sections along their group axis, ``nubar`` / ``chi``, and ``sigS`` /
    ``sig2`` along BOTH group axes — leaving ``sig0`` (a background cross section,
    NOT energy-indexed) untouched. After this the isotope obeys the canonical
    convention: group 0 = fastest, ``eg`` descending, downscatter upper-triangular.
    """
    rev = slice(None, None, -1)
    return replace(
        iso,
        eg=np.ascontiguousarray(iso.eg[rev]),
        sigC=np.ascontiguousarray(iso.sigC[:, rev]),
        sigL=np.ascontiguousarray(iso.sigL[:, rev]),
        sigF=np.ascontiguousarray(iso.sigF[:, rev]),
        sigT=np.ascontiguousarray(iso.sigT[:, rev]),
        nubar=np.ascontiguousarray(np.asarray(iso.nubar)[rev]),
        chi=np.ascontiguousarray(np.asarray(iso.chi)[rev]),
        sigS=[[_reverse_groups_2d(s) for s in order] for order in iso.sigS],
        sig2=[_reverse_groups_2d(s) for s in iso.sig2],
    )


def _extract_chi(i_temp: int, m: np.ndarray) -> np.ndarray:
    """Extract the fission spectrum (chi) from MF=6, MT=18."""
    # Find MF=6, MT=18
    i = 0
    while i < m.shape[0]:
        if m[i, 7] == 6 and m[i, 8] == 18:
            break
        i += 1
    else:
        return np.zeros(NG)

    i += 1  # move to second record
    ig2lo = int(m[i, 3])
    nw = int(m[i, 4])
    i += 1
    a, _ = _extract_n_words(nw, i, m)

    chi = np.zeros(NG)
    for j in range(nw):
        chi[ig2lo - 1 + j] = a[j]
    total = chi.sum()
    if total > 0:
        chi /= total
    return chi


def _init_scattering(
    elastic: tuple | None, n_sig0: int, n_legendre: int
) -> list[list[csr_matrix]]:
    """Initialize the ``n_legendre`` × ``n_sig0`` scattering matrix list from elastic data.

    ``n_legendre`` is the widest order any scattering section of the isotope
    stores; the ingest invents nothing and discards nothing — the solve's
    ``scattering_order`` is the only truncation (#426).
    """
    sigS: list[list[csr_matrix]] = [
        [csr_matrix((NG, NG)) for _ in range(n_sig0)] for _ in range(n_legendre)
    ]

    if elastic is None:
        return sigS

    ifrom, ito, data = elastic

    # Zero out thermal groups for elastic (they're handled by thermal scattering)
    thermal_mask = ifrom <= _IG_THRESH

    for j_lgn in range(len(sigS)):
        for i_sig0 in range(n_sig0):
            key = (j_lgn, i_sig0)
            if key in data:
                vals = data[key].copy()
                vals[thermal_mask] = 0.0
                vals += _SPARSITY_EPSILON
                # .tocsr() yields a csr_matrix at runtime; cast past the
                # csr_array inference of the untyped tocsr() body.
                sigS[j_lgn][i_sig0] = cast(
                    csr_matrix,
                    coo_matrix(
                        (vals, (ifrom - 1, ito - 1)), shape=(NG, NG)
                    ).tocsr(),
                )
    return sigS


def _accumulate_scattering(
    sigS: list[list[csr_matrix]],
    reaction: tuple,
    n_sig0: int,
    sigma_zero_independent: bool = False,
) -> None:
    """Add a scattering reaction (inelastic or thermal) into sigS, order by order."""
    ifrom, ito, data = reaction

    for j_lgn in range(len(sigS)):
        for i_sig0 in range(n_sig0):
            # Inelastic/thermal: same for all sigma-zeros (use sig0=0 data)
            src_key = (j_lgn, 0) if sigma_zero_independent else (j_lgn, i_sig0)
            if src_key in data:
                vals = data[src_key] + _SPARSITY_EPSILON
                addition = coo_matrix(
                    (vals, (ifrom - 1, ito - 1)), shape=(NG, NG)
                ).tocsr()
                sigS[j_lgn][i_sig0] = sigS[j_lgn][i_sig0] + addition
