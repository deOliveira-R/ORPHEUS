"""
WIMS-D/IAEA (WLUP) energy group structures: 69-group and 172-group.

Source
------
IAEA-TECDOC / STI/Pub/1264, "WIMS-D Library Update" (IAEA, 2007).
  - Table 11.1  "69 ENERGY GROUP STRUCTURE"   (printed p.54  / PDF p.67)
  - Table 11.2  "172 ENERGY GROUP STRUCTURE"  (printed pp.55-56 / PDF pp.68-69)
  - Table 11.3  "CORRESPONDENCE BETWEEN THE 172 AND 69 ENERGY GROUPS" (PDF p.69)

All energies are UPPER group boundaries Emax in eV. Groups are indexed 1..N from
the highest-energy group downward, matching the WIMS convention. The arrays below
store Emax for groups 1..N; the lower boundary of the bottom group is held
separately (G69_LOWER, G172_LOWER), giving N+1 edges total.

Structure summary (from the source text):
  69-group : 14 fast (9.118 keV - 10 MeV)
             13 resonance (4 eV - 9.118 keV)
             42 thermal (0 - 4 eV)          -> bottom edge is 0.0 eV
  172-group: 45 fast (9.119 keV - 19.64 MeV)
             47 resonance (4 eV - 9.119 keV)
             80 thermal (1e-5 - 4 eV)        -> bottom edge is 1.0e-5 eV

============================ TRANSCRIPTION DEFECT ============================
In the published PDF, the thermal Emax column for 172-group groups 140..172
(and the bottom edge) is CLIPPED at the page's right margin: the mantissas
render ("4.33000E", "4.00000E", ...) but the exponent characters are absent
from the PDF content stream (mantissa x1 = 537.2 pt on a 595 pt page). This is
a defect in the source document, not an extraction error.

The 34 missing exponents (groups 140..172 + bottom edge) were RECONSTRUCTED,
not transcribed, from three mutually-consistent constraints:
  1. Strict monotonic decrease of Emax with group index (forces the exponent
     to step from E-01 -> E-02 at g157, and E-02 -> E-03 at g170).
  2. Exact agreement with the 69-group thermal boundaries (which are NOT
     clipped and carry full exponents) at every shared boundary, e.g.
     g141=4.00000E-01, g146=3.00000E-01, g156=1.00000E-01, g158=8.00000E-02,
     g169=1.00000E-02, g171=5.00000E-03.
  3. The stated thermal range "10^-5 to 4 eV" fixes the bottom edge at
     1.00000E-05 eV (also printed as an explicit final table row).
The Table 11.3 correspondence independently corroborates the reconstruction
(e.g. 69 g47 -> 172 g140 with the g140/g141 boundary at 4.0e-1 eV).
Confidence: very high. If you have an electronic WLUP 172-group library
(.WIMS / GENDF header) on hand, cross-check these 34 values against it before
production use.
=============================================================================
"""

import numpy as np

from orpheus.data.energy_grid import EnergyGrid

# --------------------------------------------------------------------------
# 69-group structure  (Table 11.1).  Emax[eV] for groups 1..69.
# --------------------------------------------------------------------------
G69_EMAX = [
    1.00000E+07, 6.06550E+06, 3.67900E+06, 2.23100E+06, 1.35300E+06, 8.21000E+05,
    5.00000E+05, 3.02500E+05, 1.83000E+05, 1.11000E+05, 6.73400E+04, 4.08500E+04,
    2.47800E+04, 1.50300E+04,                                                       # fast 1-14
    9.11800E+03, 5.53000E+03, 3.51910E+03, 2.23945E+03, 1.42510E+03, 9.06899E+02,
    3.67263E+02, 1.48729E+02, 7.55014E+01, 4.80520E+01, 2.77000E+01, 1.59680E+01,
    9.87700E+00,                                                                    # resonance 15-27
    4.00000E+00, 3.30000E+00, 2.60000E+00, 2.10000E+00, 1.50000E+00, 1.30000E+00,
    1.15000E+00, 1.12300E+00, 1.09700E+00, 1.07100E+00, 1.04500E+00, 1.02000E+00,
    9.96000E-01, 9.72000E-01, 9.50000E-01, 9.10000E-01, 8.50000E-01, 7.80000E-01,
    6.25000E-01, 5.00000E-01, 4.00000E-01, 3.50000E-01, 3.20000E-01, 3.00000E-01,
    2.80000E-01, 2.50000E-01, 2.20000E-01, 1.80000E-01, 1.40000E-01, 1.00000E-01,
    8.00000E-02, 6.70000E-02, 5.80000E-02, 5.00000E-02, 4.20000E-02, 3.50000E-02,
    3.00000E-02, 2.50000E-02, 2.00000E-02, 1.50000E-02, 1.00000E-02, 5.00000E-03,  # thermal 28-69
]
G69_LOWER = 0.0  # bottom edge of group 69 ("thermal groups from 0 to 4 eV")

# --------------------------------------------------------------------------
# 172-group structure  (Table 11.2).  Emax[eV] for groups 1..172.
# Groups 140..172 exponents reconstructed -- see TRANSCRIPTION DEFECT above.
# --------------------------------------------------------------------------
G172_EMAX = [
    # fast 1-45
    1.96403E+07, 1.73325E+07, 1.49182E+07, 1.38403E+07, 1.16183E+07, 1.00000E+07,
    8.18731E+06, 6.70320E+06, 6.06531E+06, 5.48812E+06, 4.49329E+06, 3.67879E+06,
    3.01194E+06, 2.46597E+06, 2.23130E+06, 2.01897E+06, 1.65299E+06, 1.35335E+06,
    1.22456E+06, 1.10803E+06, 1.00259E+06, 9.07180E+05, 8.20850E+05, 6.08101E+05,
    5.50232E+05, 4.97871E+05, 4.50492E+05, 4.07622E+05, 3.01974E+05, 2.73237E+05,
    2.47235E+05, 1.83156E+05, 1.22773E+05, 1.11090E+05, 8.22975E+04, 6.73795E+04,
    5.51656E+04, 4.08677E+04, 3.69786E+04, 2.92830E+04, 2.73944E+04, 2.47875E+04,
    1.66156E+04, 1.50344E+04, 1.11378E+04,
    # resonance 46-92
    9.11882E+03, 7.46586E+03, 5.53085E+03, 5.00450E+03, 3.52662E+03, 3.35463E+03,
    2.24867E+03, 2.03468E+03, 1.50733E+03, 1.43382E+03, 1.23410E+03, 1.01039E+03,
    9.14242E+02, 7.48518E+02, 6.77287E+02, 4.53999E+02, 3.71703E+02, 3.04325E+02,
    2.03995E+02, 1.48625E+02, 1.36742E+02, 9.16609E+01, 7.56736E+01, 6.79040E+01,
    5.55951E+01, 5.15780E+01, 4.82516E+01, 4.55174E+01, 4.01690E+01, 3.72665E+01,
    3.37201E+01, 3.05113E+01, 2.76077E+01, 2.49805E+01, 2.26033E+01, 1.94548E+01,
    1.59283E+01, 1.37096E+01, 1.12245E+01, 9.90555E+00, 9.18981E+00, 8.31529E+00,
    7.52398E+00, 6.16012E+00, 5.34643E+00, 5.04348E+00, 4.12925E+00,
    # thermal 93-139 (exponents present in source)
    4.00000E+00, 3.38075E+00, 3.30000E+00, 2.76792E+00, 2.72000E+00, 2.60000E+00,
    2.55000E+00, 2.36000E+00, 2.13000E+00, 2.10000E+00, 2.02000E+00, 1.93000E+00,
    1.84000E+00, 1.75500E+00, 1.67000E+00, 1.59000E+00, 1.50000E+00, 1.47500E+00,
    1.44498E+00, 1.37000E+00, 1.33750E+00, 1.30000E+00, 1.23500E+00, 1.17000E+00,
    1.15000E+00, 1.12535E+00, 1.11000E+00, 1.09700E+00, 1.07100E+00, 1.04500E+00,
    1.03500E+00, 1.02000E+00, 9.96000E-01, 9.86000E-01, 9.72000E-01, 9.50000E-01,
    9.30000E-01, 9.10000E-01, 8.60000E-01, 8.50000E-01, 7.90000E-01, 7.80000E-01,
    7.05000E-01, 6.25000E-01, 5.40000E-01, 5.00000E-01, 4.85000E-01,
    # thermal 140-172 (exponents RECONSTRUCTED)
    4.33000E-01, 4.00000E-01, 3.91000E-01, 3.50000E-01, 3.20000E-01, 3.14500E-01,
    3.00000E-01, 2.80000E-01, 2.48000E-01, 2.20000E-01, 1.89000E-01, 1.80000E-01,
    1.60000E-01, 1.40000E-01, 1.34000E-01, 1.15000E-01, 1.00000E-01, 9.50000E-02,
    8.00000E-02, 7.70000E-02, 6.70000E-02, 5.80000E-02, 5.00000E-02, 4.20000E-02,
    3.50000E-02, 3.00000E-02, 2.50000E-02, 2.00000E-02, 1.50000E-02, 1.00000E-02,
    6.90000E-03, 5.00000E-03, 3.00000E-03,
]
G172_LOWER = 1.00000E-05  # bottom edge ("thermal groups from 10^-5 to 4 eV")

# --------------------------------------------------------------------------
# Table 11.3 correspondence, expressed as a condensation map.
# Each tuple is (coarse_69_group, first_172_group, last_172_group): the inclusive
# range of fine 172-groups that collapse into the given 69-group.
# NOTE: this is a NEAREST-boundary correspondence, NOT exact boundary
# coincidence -- the two structures were defined independently. See
# boundary_mismatch_report() below. In particular, 172-groups 1..5 lie ABOVE
# the 69-group ceiling (10 MeV vs 19.64 MeV) and fold into 69-group 1.
# --------------------------------------------------------------------------
_LAST_172_OF_69 = [
      8, 11, 14, 17, 22, 25, 28, 31, 33, 35, 37, 41, 43, 45, 47, 49, 51, 54,
     57, 61, 64, 67, 71, 77, 81, 84, 92, 94, 97,101,108,113,116,117,119,120,
    121,123,124,126,127,129,131,133,135,137,140,142,143,145,146,147,148,150,
    152,155,157,159,160,161,162,163,164,165,166,167,168,170,172,
]
CONDENSE_172_TO_69 = []
_start = 1
for _g69, _last in enumerate(_LAST_172_OF_69, start=1):
    CONDENSE_172_TO_69.append((_g69, _start, _last))
    _start = _last + 1


def edges(emax, lower):
    """Return the N+1 monotone-decreasing edge array for a structure."""
    return list(emax) + [lower]


# Canonical descending EnergyGrid instances (the energy-condensation targets).
# The WIMS arrays store group 1 = fastest (the WIMS convention); read as an
# edges() array they are ALREADY in the canonical ORPHEUS descending orientation
# (eg[0] = highest E, group 0 = fastest), so no reversal is needed -- only the
# 1-based WIMS label maps to the 0-based ORPHEUS group index.  CONDENSE_172_TO_69
# (above) stays the derivation-validation oracle for fine.overlap_to(coarse).
WIMS_69 = EnergyGrid(np.asarray(edges(G69_EMAX, G69_LOWER)))
WIMS_172 = EnergyGrid(np.asarray(edges(G172_EMAX, G172_LOWER)))


def validate():
    """Cheap structural V&V. Raises AssertionError on any inconsistency."""
    assert len(G69_EMAX) == 69
    assert len(G172_EMAX) == 172
    assert len(CONDENSE_172_TO_69) == 69
    for name, e, lo in (("69", G69_EMAX, G69_LOWER), ("172", G172_EMAX, G172_LOWER)):
        seq = edges(e, lo)
        assert all(seq[i] > seq[i + 1] for i in range(len(seq) - 1)), \
            f"{name}-group edges not strictly decreasing"
    flat = [g for (_, a, b) in CONDENSE_172_TO_69 for g in range(a, b + 1)]
    assert flat == list(range(1, 173)), "condensation map is not a clean partition of 1..172"
    return True


def boundary_mismatch_report(rel_tol=5e-4):
    """
    For each 69-group, compare its UPPER edge to the upper edge of the first
    172-group assigned to it. Returns list of (g69, first172, e69, e172, reldiff)
    for boundaries differing by more than rel_tol. Demonstrates that condensation
    via Table 11.3 is approximate at these boundaries.
    """
    out = []
    for g69, a, _b in CONDENSE_172_TO_69:
        e69 = G69_EMAX[g69 - 1]
        e172 = G172_EMAX[a - 1]
        rel = abs(e69 - e172) / e69
        if rel > rel_tol:
            out.append((g69, a, e69, e172, rel))
    return out


if __name__ == "__main__":
    validate()
    print("validate(): OK")
    print(f"69-group: {len(G69_EMAX)} groups, ceiling {G69_EMAX[0]:.4E} eV, floor {G69_LOWER:.4E} eV")
    print(f"172-group: {len(G172_EMAX)} groups, ceiling {G172_EMAX[0]:.4E} eV, floor {G172_LOWER:.4E} eV")
    mism = boundary_mismatch_report()
    print(f"\n{len(mism)} coarse boundaries differ from the nearest fine boundary by >0.05%:")
    for g69, a, e69, e172, rel in mism:
        print(f"  69 g{g69:<2d} {e69:.5E} eV  vs  172 g{a:<3d} {e172:.5E} eV  ({rel*100:.3f}%)")
