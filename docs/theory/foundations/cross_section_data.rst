.. _theory-cross-section-data:

==============================================
Cross-Section Data Pipeline
==============================================

.. contents:: Contents
   :local:
   :depth: 3


Key Facts
=========

**Read this before modifying the cross-section pipeline.**

- 421-group microscopic XS from GENDF (GXS) files → HDF5 via ``data/micro_xs/``
- **Canonical energy-group convention**: group index ``0`` = **fastest**
  (highest energy), ``eg`` strictly **descending**, ``SigS[g_from, g_to]``
  downscatter in the **upper triangle**. Enforced once at the GENDF ingest
  (NJOY is the opposite); see :ref:`canonical-group-convention`.
- ``Isotope`` dataclass: sig_t, sig_c, sig_f, sig_el, sig_inel, nu, chi (421 groups)
- Sigma-zero iteration: ``orpheus/data/macro_xs/sigma_zeros.py`` — self-shielding
- ``Mixture`` dataclass: macroscopic XS with ``SigS[l][g_from, g_to]`` convention
- Consistency: :math:`\Sigma_t = \Sigma_c + \Sigma_f + \sum_{g'} \Sigma_{s,g \to g'}`
- ``load_isotope()`` auto-selects HDF5 or fallback .m parser
- Verification uses synthetic XS from ``orpheus/derivations/common/xs_library.py`` (regions A/B/C/D), NOT this pipeline


Overview
========

Every solver in ORPHEUS relies on multi-group microscopic cross sections
for the 12 nuclides in the 421-energy-group JEFF-3.1 library.  This
chapter documents the complete data pipeline from the authoritative IAEA
source files to the internal ``Isotope`` dataclass:

1. **GENDF format** — the IAEA distributes processed nuclear data in
   the GENDF (Groupwise ENDF) format, which uses fixed-width 80-column
   records inherited from punched-card conventions.
2. **Parsing** — the ``gendf.py`` module reads GENDF files directly,
   bypassing the MATLAB CSV intermediary.
3. **Scattering assembly** — elastic, inelastic, and thermal scattering
   matrices are combined with careful treatment of thermal-group
   boundaries.
4. **HDF5 serialisation** — the parsed data is stored in compressed HDF5
   files for fast loading at runtime.
5. **Loading** — ``load_isotope()`` provides a uniform API that
   auto-selects the HDF5 backend when available.

The 12 nuclides in the library are:

.. list-table::
   :header-rows: 1
   :widths: 15 15 15 20

   * - Nuclide
     - GXS File
     - Temperatures (K)
     - Sigma-zeros
   * - H-1
     - ``H_001.GXS``
     - 294, 350, 400, 450, 500, 550, 600, 650
     - 1
   * - O-16
     - ``O_016.GXS``
     - 294, 600, 900, 1200, 1500, 1800
     - 6
   * - B-10
     - ``B_010.GXS``
     - 294, 600, 900, 1200
     - 4
   * - B-11
     - ``B_011.GXS``
     - 294, 600, 900, 1200
     - 4
   * - Na-23
     - ``NA023.GXS``
     - 294, 600, 900, 1200
     - 4
   * - U-235
     - ``U_235.GXS``
     - 294, 600, 900, 1200, 1500, 1800
     - 10
   * - U-238
     - ``U_238.GXS``
     - 294, 600, 900, 1200, 1500, 1800
     - 10
   * - Zr-90
     - ``ZR090.GXS``
     - 294, 600, 900, 1200
     - 4
   * - Zr-91
     - ``ZR091.GXS``
     - 294, 600, 900, 1200
     - 4
   * - Zr-92
     - ``ZR092.GXS``
     - 294, 600, 900, 1200
     - 4
   * - Zr-94
     - ``ZR094.GXS``
     - 294, 600, 900, 1200
     - 4
   * - Zr-96
     - ``ZR096.GXS``
     - 294, 600, 900, 1200
     - 4


.. _canonical-group-convention:

The Canonical Energy-Group Convention
=====================================

.. note::

   **This section is the single source of truth for energy-group
   ordering in ORPHEUS.** Every page, solver, and data array obeys it;
   where another page restates it, that page defers here.

ORPHEUS uses **one** energy-group convention everywhere. State it once,
precisely:

- **Group index** :math:`g = 0` is the **fastest** group (highest
  energy); :math:`g = G-1` is the **slowest** (thermal). Equivalently,
  groups are ordered fast → thermal as the index increases.
- **The group-boundary array** :attr:`Isotope.eg <orpheus.data.micro_xs.isotope.Isotope>`
  (length :math:`G+1`) is strictly **descending**:
  :math:`\mathrm{eg}[0] > \mathrm{eg}[1] > \cdots > \mathrm{eg}[G]`, so
  :math:`\mathrm{eg}[0]` is the highest energy bound and
  :math:`\mathrm{eg}[G]` the lowest. Group :math:`g` spans the interval
  :math:`[\mathrm{eg}[g+1],\, \mathrm{eg}[g]]`.
- **The scattering matrix** :math:`\Sigma_{\mathrm{s},\ell}[g_{\text{from}},
  g_{\text{to}}]` therefore carries **downscatter** (a neutron losing
  energy: :math:`g_{\text{to}} > g_{\text{from}}`) in the **upper
  triangle**. **Upscatter** (energy gain: :math:`g_{\text{to}} <
  g_{\text{from}}`, the lower triangle) is physically possible only among
  the **thermal groups**, where neutrons are in near-equilibrium with the
  thermal motion of the medium.

.. important::

   This is the ``[g_from, g_to]`` row-source convention shared with the
   rest of the codebase: ``SigS[l]`` stores
   :math:`\Sigma_{\mathrm{s},\ell}[g_{\text{from}}, g_{\text{to}}]`, so
   the in-scatter source contracts the **transpose**,
   :math:`q_{\mathrm{s}} = \Sigma_{\mathrm{s}}^{\mathsf T}\phi`. This
   energy-index convention is **orthogonal** to the *array storage
   layout* :math:`(N, n_g, n_x, n_y)` (the *axis ordering* within each
   array, documented in :ref:`the SN index-convention page
   <theory-sn-index-convention>`) — this section fixes *which physical
   group each index labels*, not *which array axis the group sits on*.

Enforcement: one normalisation at the data-ingest boundary
----------------------------------------------------------

The convention is enforced **exactly once**, at the boundary where the
foreign NJOY/GENDF data enters ORPHEUS. The IAEA GENDF files store the
**opposite** order — group 0 = thermal, energies **ascending** — a
historical NJOY artefact. The single normalisation
:func:`_to_canonical_group_order <orpheus.data.micro_xs.gendf>`
(with :func:`_reverse_groups_2d <orpheus.data.micro_xs.gendf>` for the
two-axis :math:`[g_{\text{from}}, g_{\text{to}}]` matrix channels)
reverses every group-indexed array of an :class:`~orpheus.data.micro_xs.isotope.Isotope`
**once on ingest**, inside :func:`convert_gxs <orpheus.data.micro_xs.gendf>`:

- the vector cross sections (:math:`\sigma_t, \sigma_c, \sigma_f,
  \sigma_L`) are reversed along their group axis;
- :math:`\bar\nu` and :math:`\chi` are reversed;
- the scattering matrices ``sigS`` and the :math:`(n,2n)` matrix
  ``sig2`` are reversed along **both** group axes (moving downscatter
  from the lower to the upper triangle);
- :math:`\sigma_0` — a *background* cross section, **not** energy-indexed
  — is left untouched;
- ``eg`` is reversed to descending.

After this single flip, **every downstream consumer is order-transparent**:
:class:`~orpheus.data.macro_xs.mixture.Mixture`,
:func:`~orpheus.data.macro_xs.mixture.compute_macro_xs`, and every solver
read the canonical fast-first arrays directly and never re-order. The
``.h5`` caches are written **after** the flip, so cached data is
fast-first as well.

.. warning::

   The thermal-cutoff constant ``_IG_THRESH = 95`` in
   ``orpheus/data/micro_xs/gendf.py`` indexes the **native NJOY
   (ascending)** order, because it is applied during *extraction*
   (:func:`_init_scattering <orpheus.data.micro_xs.gendf>`,
   ``thermal_mask = ifrom <= _IG_THRESH``), **before**
   :func:`_to_canonical_group_order` runs. In the *stored* canonical
   arrays the same thermal boundary lives at a **high** index (near
   :math:`G - 95 \approx 326`), not at index 95. See
   :ref:`thermal-group-boundary` for the full extraction-vs-storage
   distinction.

Why the fast-first convention (rationale)
-----------------------------------------

Three reasons, in order of weight:

1. **It converges the codebase onto one convention.** Before the flip
   the GENDF pipeline was the *lone* group-0-thermal outlier. The
   synthetic verification cross sections
   (:mod:`~orpheus.derivations.common.xs_library`), the Sood benchmark
   registry (:doc:`/theory/references/sood_registry`), the :mod:`~orpheus.diffusion`
   solver, the test fixtures, and **every theory page** were already
   group-0-fast. Reversing the GENDF data on ingest makes the production
   nuclear data agree with all of them for the first time, removing a
   per-source "which way does this array run?" hazard.

2. **It is physics-identical — the multigroup solve is order-agnostic.**
   ORPHEUS does **not** run an energy-group Gauss–Seidel sweep (a
   downscatter cascade that *would* depend on the spectral ordering — and
   which does not exist in any ORPHEUS solver; see
   :ref:`group-vs-octant-terminology`). The in-scatter source is the
   **full contraction** over all source groups,

   .. math::
      :label: in-scatter-full-contraction

      q_{\mathrm{s},g} = \sum_{g'} \Sigma_{\mathrm{s},\, g' \to g}\, \phi_{g'},

   which is invariant under any permutation of the group index applied
   consistently to :math:`\Sigma_{\mathrm{s}}`, :math:`\phi`, and
   :math:`\chi`. A permutation-invariance gate confirms this: under a
   group reversal the eigenvalue :math:`\keff` is invariant and the flux
   simply reverses with the index. The flip is therefore a relabelling,
   not a change of the solved problem.

3. **Fast-first gives the natural downscatter sweep direction.** With
   :math:`g` increasing from fast to thermal, slowing-down flows in the
   direction of increasing index, so the downscatter source matrix is
   upper-triangular and the spectrum is read top-to-bottom fast → thermal
   — the way reactor-physics texts present the slowing-down equation.

.. vv-status: in-scatter-full-contraction documented
.. (vv-status rationale) representational: this is the group-coupling
.. definition that makes the order-agnosticism argument, not a solver
.. claim. The order-invariance it asserts is pinned by the permutation
.. gate referenced above (keff invariant, flux reversed); the equation
.. itself transcribes the standard multigroup in-scatter source.


.. _group-vs-octant-terminology:

Terminology: "group" — energy group vs. octant group
====================================================

The word **"group"** is overloaded in transport, and the two meanings
live on **orthogonal axes**. The canonical convention above is entirely
about the *energy* axis; do **not** conflate it with the *angular* sense
the SN :term:`sweep` uses.

.. list-table::
   :header-rows: 1
   :widths: 18 41 41

   * - Axis
     - "group" means …
     - Where it appears
   * - **Energy**
     - One of the :math:`G` energy bins of the multigroup spectrum
       (:ref:`canonical-group-convention`). Index 0 = fastest.
     - ``eg``, ``SigS[g_from, g_to]``, :math:`\phi_g`, every cross
       section. **This** is what the canonical convention orders.
   * - **Angular (SN sweep)**
     - An **octant group**: a set of :term:`ordinate` directions that share a
       sweep direction (a wavefront sweeps the spatial mesh once per
       octant group).
     - The ``inner_schedule="gauss_seidel"`` sweep schedule
       (:class:`SweepSchedule <orpheus.sn.loss_representation.sweep_schedule.SweepSchedule>`).

In the SN solver, selecting ``inner_schedule="gauss_seidel"`` (see
:func:`~orpheus.sn.solver._select_si_splitting` and the reified
splitting matrix
:class:`~orpheus.sn.operators.scheduled_invertible.ScheduledInvertibleOperator`)
builds an **octant-group / boundary** Gauss–Seidel: it folds the
reflective boundary operator :math:`B` **into** a multi-dimensional
wavefront sweep,
re-reflecting each octant group's outgoing reflective faces between octant
sweeps so a later octant reads the fresh current-iterate inflow. This is a
purely **angular** acceleration of the within-group source iteration; its
converged **bulk** is identical to the plain (Jacobi) sweep — only the
iteration rate changes (it is splitting-invariant in the sense of
``vv-principles`` Mode 9). It does **not** touch the energy index, and it
does **not** cascade across the spectrum.

.. note::

   ⛔ *"Converged fixed point is identical"* until 2026-08-15 (#344).
   The **bulk** claim is right and is what every gate measures; the
   **trace** claim is not, whenever :math:`A = L+C-S-B` is singular —
   which it exactly is on a diamond-difference mesh closing
   :math:`\ge 2` reflective axis pairs.  There the two schedules
   legitimately return different **members** of a solution manifold
   (``[M]`` differing by :math:`0.124184`, **100 %** of it inside
   :math:`\ker A`), and the solver gauges the returned trace onto the
   canonical member.  See :ref:`sn-loss-kernel-gauge`.

.. note::

   **An *energy-group* Gauss–Seidel — a downscatter cascade that solves
   group** :math:`g` **using the already-updated fluxes of the faster
   groups** :math:`g' < g` **— does NOT exist in ORPHEUS today.** Every
   solver treats the in-scatter source as the full contraction
   :eq:`in-scatter-full-contraction` over *all* source groups (lagged
   from the previous outer iterate), with **no** sweep ordering imposed
   on the energy axis. This is exactly *why* the energy-group ordering is
   free to be a pure relabelling (rationale point 2 above): there is no
   energy-sweep whose direction the convention could change. A reader
   must therefore never read the octant ``gauss_seidel`` schedule as a
   spectral sweep — the two are unrelated.


The GENDF Format
=================

Source
------

The GENDF (Groupwise ENDF) files are obtained from the IAEA Nuclear Data
Services:

   https://www-nds.iaea.org/ads/adsgendf.html

These are the 421-group JEFF-3.1 processed nuclear data files.  Each file
contains all reaction cross sections and transfer matrices for one
nuclide at multiple temperatures.


Record Layout
--------------

Every line in a GENDF file is exactly 80 characters wide, following the
ENDF-6 format :cite:`ENDF102` inherited from the punched-card era:

.. code-block:: text

   Columns  1-66: 6 data fields, 11 characters each
   Columns 67-70: MAT number (material identifier)
   Columns 71-72: MF  number (file type)
   Columns 73-75: MT  number (reaction type)
   Columns 76-80: line sequence number

For example, the first data record of H-1 looks like:

.. code-block:: text

    1.001000+3 9.991673-1          0          1         -1          1 125 1451    1
   |----11----|----11----|----11----|----11----|----11----|----11----| MAT|MF|MT |SEQ|


Compact Float Notation
-----------------------

Data fields use a compact Fortran notation where the ``E`` in scientific
notation is omitted.  The exponent sign immediately follows the mantissa:

.. code-block:: text

   1.001000+3  →  1.001000E+3  =  1001.0
   9.991673-1  →  9.991673E-1  =  0.9991673
   2.407191-7  →  2.407191E-7  =  2.407191×10⁻⁷

The parser ``_parse_gendf_field`` in ``gendf.py`` handles this
by inserting ``E`` before any ``+`` or ``-`` sign that follows a digit:

.. code-block:: python

   s = re.sub(r"(\d)([+-])", r"\1E\2", s)
   return float(s)


MF and MT Numbers
------------------

The MF (file) and MT (reaction) numbers identify the type of data in
each section:

**MF=1 — General information:**

.. list-table::
   :header-rows: 1
   :widths: 10 40

   * - MT
     - Content
   * - 451
     - Header: temperatures, sigma-zero base points, energy group
       boundaries (422 values for 421 groups)

**MF=3 — Cross sections** (sigma-zero dependent, one value per group):

.. list-table::
   :header-rows: 1
   :widths: 10 40

   * - MT
     - Reaction
   * - 1
     - Total cross section (does **not** include upscattering)
   * - 18
     - Fission
   * - 102
     - Radiative capture :math:`(n,\gamma)`
   * - 107
     - :math:`(n,\alpha)`
   * - 452
     - Total :math:`\bar{\nu}` (average neutrons per fission)

**MF=6 — Transfer matrices** (group-to-group scattering):

.. list-table::
   :header-rows: 1
   :widths: 10 40

   * - MT
     - Reaction
   * - 2
     - Elastic scattering (sigma-zero dependent)
   * - 16
     - :math:`(n,2n)` reaction
   * - 18
     - Fission spectrum :math:`\chi(g)`
   * - 51–91
     - Discrete inelastic scattering levels
   * - 221
     - Free-gas thermal scattering
   * - 222
     - Thermal scattering for H bound in water (:math:`S(\alpha,\beta)`)


MF=3 Record Structure
----------------------

Each MF=3 section begins with a header record followed by per-group
data records.  The structure for a section with :math:`N_\ell` Legendre
components and :math:`N_{\sigma_0}` sigma-zero values is:

.. code-block:: text

   Record 1 (section header):
     [ZA, AWR, NL, N_sig0, LRFLAG, NG, MAT, MF, MT, 1]

   For each group g = 1, ..., NG:
     Record (group header):
       [TEMP, 0, NL, N_sig0, NW, IG, MAT, MF, MT, line]
     Record(s) (data):
       NW = 2 × NL × N_sig0 words packed 6 per line

The first half of the NW words contains flux weights; the second half
contains the cross-section values organised as:

.. math::

   a[N_\ell N_{\sigma_0} + 1 : N_\ell N_{\sigma_0} + N_{\sigma_0}]
   = \sigma_{x,g}(\sigma_{0,1}), \ldots, \sigma_{x,g}(\sigma_{0,N_{\sigma_0}})

This is the Legendre-0 component for each sigma-zero.  Higher Legendre
components follow in the same block.


MF=6 Record Structure
----------------------

Transfer matrices in MF=6 are stored per source group in a sparse
representation.  For each source group :math:`g`:

.. code-block:: text

   Record (group header):
     [TEMP(?), 0, NG2, IG2LO, NW, IG, MAT, MF, MT, line]
   Record(s) (data):
     NW words packed 6 per line

where:

- ``NG2`` — number of secondary (target) groups with non-zero values
- ``IG2LO`` — 1-based index of the lowest non-zero target group
- ``NW`` — total words to read (includes flux weights)
- ``IG`` — 1-based source group index

The data layout per source group is:

1. **Flux weights**: :math:`N_\ell \times N_{\sigma_0}` values (skipped)
2. **Transfer values**: for each target group from ``IG2LO`` to
   ``IG2LO + NG2 - 2``, and for each sigma-zero and Legendre order:

   .. code-block:: text

      for i_to = IG2LO to IG2LO + NG2 - 2:
          for i_sig0 = 1 to N_sig0:
              for i_lgn = 1 to N_lgn:
                  sigma_s(IG → i_to, Legendre=i_lgn, sig0=i_sig0)


Scattering Matrix Assembly
===========================

The scattering matrix :math:`\Sigma_{\mathrm{s},\ell}^{(\sigma_0)}`
is assembled from three separate GENDF sections.  This is one of the
most delicate parts of the pipeline.


.. _thermal-group-boundary:

Thermal-Group Boundary
-----------------------

The free-atom elastic scattering model breaks down below
:math:`E \approx 4` eV, where the target atoms are bound in a lattice or
molecule (thermal motion affects scattering), so the assembly zeroes the
elastic kernel in the thermal range and replaces it with a bound-atom
thermal kernel.

.. important::

   **This boundary is expressed in the NATIVE NJOY group index, not the
   ORPHEUS canonical one.** All of the scattering assembly below runs on
   the raw GENDF arrays *before* the ingest flip
   (:func:`_to_canonical_group_order <orpheus.data.micro_xs.gendf>`,
   :ref:`canonical-group-convention`), where group 0 = thermal and the
   index increases with energy. In that native order the thermal cutoff
   is the constant ``_IG_THRESH = 95``: native groups :math:`g \le 95`
   (:math:`E \lesssim 4` eV) are the thermal range. After the flip, the
   thermal groups occupy the **high-index** end of the canonical arrays
   (near :math:`G - 95 \approx 326`), consistent with group 0 = fastest.
   Read every ``g \le 95`` / ``g > 95`` reference in this section as
   *native-order* indices.

The GENDF file provides two models:

- **MT=2** — free-atom elastic scattering (valid above ~4 eV)
- **MT=221** — free-gas thermal scattering (all isotopes except H-1)
- **MT=222** — :math:`S(\alpha,\beta)` thermal scattering for
  H bound in water (H-1 only)


Assembly Algorithm
-------------------

The scattering matrix is built in four stages:

1. **Elastic (MT=2)**: Extract the elastic scattering transfer matrix.
   **Zero out** all entries where the source group :math:`g \le 95` (the
   thermal range, **native NJOY index**), because thermal scattering
   replaces elastic in that range.  Add :math:`10^{-30}` to all values
   (matching the MATLAB convention to avoid exact zeros in sparse
   matrices).

   .. code-block:: python

      vals[thermal_mask] = 0.0
      vals += 1e-30
      sigS[lgn][sig0] = sparse(ifrom-1, ito-1, vals, NG, NG)

2. **Inelastic (MT=51–91)**: For each discrete inelastic level that
   exists, extract the transfer matrix and **add** it to sigS.
   Inelastic scattering is sigma-zero independent (same values for all
   sigma-zero variants), so the first sigma-zero's data is used for all.

3. **Thermal (MT=221 or MT=222)**: Extract the thermal scattering
   kernel and **add** it to sigS.  This replaces the zeroed elastic
   entries in the thermal range.  Like inelastic, thermal scattering
   is sigma-zero independent.

4. The final scattering matrix structure is a list of lists:
   ``sigS[legendre_order][sig0_index]``, each a ``csr_matrix(NG, NG)``.
   Three Legendre orders (P0, P1, P2) are always stored.

.. important::

   The elastic scattering data for groups :math:`g > 95` (epithermal
   and fast, **native NJOY index**) **does** depend on sigma-zero.  Each
   sigma-zero variant has different elastic values at these groups.  For
   groups :math:`g \le 95`, the elastic is zeroed and replaced by the
   sigma-zero-independent thermal kernel.


Reactions Not Included: :math:`(n,2n)`, :math:`(n,3n)`, :math:`(n,4n)`
-----------------------------------------------------------------------

The GENDF files for heavy isotopes (U-235, U-238, Pu-239, ...) contain
MF=6 scattering entries for three multiplicity-changing reactions that
ORPHEUS **does not extract** into the scattering matrix:

.. csv-table::
   :header: MT, Reaction, Threshold, ENDF name
   :widths: 8, 15, 15, 62

   16, ":math:`(n,2n)`", ~6–8 MeV, neutron-induced two-neutron emission
   17, ":math:`(n,3n)`", ~11–14 MeV, neutron-induced three-neutron emission
   37, ":math:`(n,4n)`", ~20 MeV+, neutron-induced four-neutron emission

The current scattering-matrix assembly loop at
``orpheus/data/micro_xs/gendf.py:301`` only iterates over
MT=51..91 (discrete inelastic levels plus continuum inelastic at
MT=91), matching the original MATLAB ``convertCSVtoM.m`` that ORPHEUS
was ported from.

Why this is (currently) acceptable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Two reasons, physical and pragmatic:

1. **Threshold energies are far above the thermal and epithermal
   regime.** A thermal-spectrum reactor has almost no flux above
   6 MeV (the fission spectrum rolls off as :math:`\chi(E) \propto
   e^{-E/a}\sinh\sqrt{bE}`, effectively dead by 10 MeV). The rate
   density for :math:`(n,xn)` at x≥2 is
   :math:`\int_{E_{\mathrm{th}}}^\infty \phi(E)\,\sigma_{(n,xn)}(E)\,dE`
   and both the flux and the cross section are negligible in the
   integration window. Quantitatively, for a PWR-like spectrum, the
   :math:`(n,2n)` rate on U-238 is below :math:`10^{-6}` of the
   absorption rate — below the truncation noise of the 421-group
   multi-group flux itself.

2. **Neutron-multiplication accounting would need to change
   consistently.** Including :math:`(n,xn)` correctly requires more
   than adding an MF=6 block to the scattering matrix. Because these
   reactions change neutron multiplicity, they must be accounted for
   separately from fission in the balance equation (they are *not*
   fission, so they do not carry :math:`\chi` or :math:`\nu`, but
   they *do* produce excess neutrons). ORPHEUS's current balance
   equation assumes a 1-in-1-out scattering model; retrofitting
   :math:`(n,xn)` cleanly means either

   - treating them as sources with a separate multiplicity factor
     (the "``nu_n_xn``" convention in MCNP/Serpent/OpenMC), or
   - folding them into the scattering matrix with an effective
     :math:`\Sigma_{\mathrm{s}}^{(n,xn)}` that scales by the
     multiplicity — which breaks the :math:`\Sigma_t`
     consistency relation :math:`\Sigma_t = \Sigma_c + \Sigma_f +
     \sum_{g'} \Sigma_{s,g\to g'}` unless :math:`\Sigma_c` is
     simultaneously adjusted.

   Doing either correctly is a data-pipeline-wide change, not a
   localized extension.

When this exclusion would matter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Three application regimes **do** need :math:`(n,xn)`:

- **Fast reactors** (SFR, LFR, GFR): fission spectrum is harder,
  5–15 % of flux above 1 MeV, :math:`(n,2n)` on U-238 and Pu-240
  contributes measurably to the fast-group balance.
- **Fusion blankets**: 14 MeV D-T source neutrons sit directly in
  the peak of the :math:`(n,2n)`/:math:`(n,3n)` cross section for
  Li, Be, and Pb — these reactions are the *whole point* of a
  breeding blanket.
- **High-energy shielding / accelerator-driven systems**: spallation
  neutron sources produce a significant population above 20 MeV,
  where :math:`(n,4n)` on heavy targets is non-negligible.

None of these are current ORPHEUS use cases. The V&V suite
(:doc:`/theory/verification/index`) exclusively verifies thermal-spectrum
analytical benchmarks; synthetic cross sections in
``orpheus/derivations/common/xs_library.py`` do not include an
:math:`(n,xn)` term.

Implementation sketch (deferred)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If / when ORPHEUS expands to fast or fusion applications, the
retrofit would touch:

1. **``orpheus/data/micro_xs/gendf.py``** — extend the scattering
   loop at line 281 to also iterate MT in ``(16, 17, 37)`` and stash
   the extracted block with a ``multiplicity`` attribute (``2``,
   ``3``, ``4`` respectively).

2. **``orpheus/data/micro_xs/isotope.py``** — add a
   ``sig_n_xn: dict[int, np.ndarray]`` field to ``Isotope``.

3. **``orpheus/data/macro_xs/mixture.py``** — decide whether
   :math:`(n,xn)` enters the ``SigS`` matrix (with a multiplicity
   factor baked in, losing the :math:`\Sigma_t` check) or lives as
   an explicit source term in the balance equation.

4. **Every transport solver** (``cp``, ``sn``, ``moc``, ``mc``,
   ``diffusion``) — the balance residual at the cell level must
   account for the multiplicity.

5. **V&V** — add an L2 benchmark against a published fast-reactor
   eigenvalue (e.g., the GODIVA or Jezebel ICSBEP criticals) where
   :math:`(n,xn)` measurably shifts :math:`k_{\mathrm{eff}}`.

Tracked in GitHub issue `#63
<https://github.com/deOliveira-R/ORPHEUS/issues/63>`_ as a
"status: impl, intentional exclusion, documentation complete"
line item.


Total Cross Section
====================

The total cross section :math:`\Sigma_{\mathrm{t},g}(\sigma_0)` is
**computed** from the components rather than read from MF=3 MT=1:

.. math::
   :label: sigT-computed

   \sigma_{\mathrm{t},g}(\sigma_0)
   = \sigma_{\mathrm{c},g}(\sigma_0)
   + \sigma_{\mathrm{f},g}(\sigma_0)
   + \sigma_{\alpha,g}(\sigma_0)
   + \sum_{g'} \sigma_{\mathrm{s},0,g \to g'}(\sigma_0)
   + \sum_{g'} \sigma_{\mathrm{2n},g \to g'}

This approach is used because:

1. MF=3 MT=1 does **not** include upscattering (stated in the MATLAB
   source: *"note that mf=3 mt=1 does not include upscatters"*).
2. Computing from components ensures self-consistency between the total
   and the reaction rates used by the solver.


.. _sigT-consistency:

sigT Consistency Issue (Historical)
-------------------------------------

.. warning::

   The legacy MATLAB ``.m`` data files contain a systematic discrepancy
   in the stored ``sigT`` values.  This section documents the issue for
   future reference.

The MATLAB ``convertCSVtoM.m`` script computes ``sigT`` from
full-precision intermediate variables and writes it with ``%13.6e``
format (6 decimal places in scientific notation).  It independently
truncates all component cross sections (sigC, sigF, sigS) to the same
format.

When the ``.m`` file is loaded and ``sigT`` is **recomputed** from the
stored (truncated) components, the result differs from the stored
``sigT`` by a **constant offset** of 10–30 barns for heavy isotopes:

.. list-table::
   :header-rows: 1
   :widths: 15 20 20 15

   * - Isotope
     - .m sigT[0,0]
     - Recomputed
     - Offset
   * - U-235 (294K)
     - 15,523.0
     - 15,504.2
     - 18.8
   * - U-238 (600K)
     - 108.14
     - 77.87
     - 30.27
   * - Zr-90 (600K)
     - (offset)
     - (recomputed)
     - 8.3
   * - H-1 (294K)
     - (matches)
     - (matches)
     - ~0

The offset is **constant across all energy groups and all sigma-zero
rows** for a given isotope/temperature.  Isotopes with :math:`N_{\sigma_0} = 1`
(like H-1) show no discrepancy.

**Root cause**: The full-precision ``sigS`` row sums differ from the
truncated-then-resummed version.  Although each individual truncation
error is :math:`O(10^{-7})` relative, the scattering matrices have
thousands of non-zero entries per row at resonance energies, and the
accumulation of truncation errors is significant.

**Impact on sigma-zero iterations**: The sigma-zero solver interpolates
``sigT(\sigma_0)`` from the tabulated values.  Using the GENDF-computed
``sigT`` instead of the ``.m`` file's stored ``sigT`` shifts the
converged sigma-zeros, which propagates to different interpolated
cross sections and ultimately a ~0.4% shift in the PWR-like mixture's
:math:`\kinf` (1.01771 vs 1.01357).


HDF5 Storage Format
=====================

Each element is stored in a single HDF5 file (e.g., ``U_235.h5``) with
one group per temperature:

.. code-block:: text

   /{temp_K}K/
       @aw          : scalar (atomic weight in amu)
       @temp        : scalar (temperature in K)
       eg           : (NG+1,)    — energy group boundaries (eV)
       sig0         : (N_sig0,)  — sigma-zero base points (barns)
       sigC         : (N_sig0, NG) — radiative capture
       sigL         : (N_sig0, NG) — (n,alpha)
       sigF         : (N_sig0, NG) — fission
       sigT         : (N_sig0, NG) — total
       nubar        : (NG,) — average neutrons per fission
       chi          : (NG,) — fission spectrum (normalised to 1)
       sig2/
           row      : (nnz,) int32  — COO row indices
           col      : (nnz,) int32  — COO column indices
           data     : (nnz,) float64 — COO values
       sigS/
           L{j}_S{k}/          — Legendre order j, sigma-zero k
               row  : (nnz,)
               col  : (nnz,)
               data : (nnz,)

Dense arrays use gzip compression (level 4).  Sparse matrices are stored
as COO triplets to avoid scipy-specific formats.


File Sizes
-----------

.. list-table::
   :header-rows: 1
   :widths: 15 15 15

   * - Element
     - Temperatures
     - HDF5 Size (MB)
   * - H-1
     - 8
     - 12.3
   * - U-235
     - 6
     - 50.0
   * - U-238
     - 6
     - 37.8
   * - O-16
     - 6
     - 10.8
   * - Zr isotopes (×5)
     - 4 each
     - ~11 each


Data Loading API
=================

The ``load_isotope`` function provides a uniform API:

.. code-block:: python

   from orpheus.data.micro_xs import load_isotope

   iso = load_isotope("U_235", 600)
   # iso.sigC — shape (10, 421), capture XS for 10 sigma-zeros
   # iso.sigS[0][0] — csr_matrix(421, 421), P0 scattering at sig0=0
   # iso.eg — shape (422,), energy group boundaries in eV

The loader reads from the HDF5 files in ``data/micro_xs/{name}.h5``.


Conversion Script
------------------

To regenerate the HDF5 files from the GENDF sources:

.. code-block:: bash

   cd orpheus/data/micro_xs
   python convert_gxs_to_hdf5.py

This processes all 12 ``.GXS`` files and writes the corresponding
``.h5`` files.  Runtime is approximately 2–3 minutes on a modern
laptop.


Validation
===========

The HDF5 data pipeline was validated by running both homogeneous
reactor cases and comparing against the MATLAB reference:

.. list-table::
   :header-rows: 1
   :widths: 30 20 20 10

   * - Case
     - Python :math:`\kinf`
     - MATLAB :math:`\kinf`
     - Match
   * - Aqueous (H₂O + U-235, 294K)
     - 1.03596
     - 1.03596
     - Yes
   * - PWR-like (UO₂ + Zry + H₂O+B, 600K)
     - 1.01357
     - 1.01357
     - Yes

Per-component validation for H-1 at 294K (1 sigma-zero, simplest case):

.. list-table::
   :header-rows: 1
   :widths: 20 20 15

   * - Quantity
     - Max diff (GXS vs .m)
     - Status
   * - ``aw``
     - 0
     - Exact
   * - ``eg``
     - 0
     - Exact
   * - ``sigC``
     - 0
     - Exact
   * - ``sigF``
     - 0
     - Exact
   * - ``nubar``
     - 0
     - Exact
   * - ``chi``
     - 0
     - Exact
   * - ``sigS[0][0]`` (row sums)
     - 0
     - Exact
   * - ``sigS[0][0]`` (nnz)
     - 77,627 = 77,627
     - Exact

Per-component validation for U-235 at 294K (10 sigma-zeros):

.. list-table::
   :header-rows: 1
   :widths: 20 20 15

   * - Quantity
     - Max diff (GXS vs .m)
     - Status
   * - ``sigC``
     - 0
     - Exact
   * - ``sigF``
     - 0
     - Exact
   * - ``nubar``
     - 0
     - Exact
   * - ``sigS[0][0]`` (row sums)
     - :math:`9.6 \times 10^{-7}`
     - Negligible
   * - ``sig2`` (nnz)
     - 6,067 = 6,067
     - Exact


.. _emission-spectrum-simplex-law:

The Fission Emission Spectrum as a Validated Value-Object
==========================================================

The fission emission spectrum :math:`\chi_g` — the probability that a
neutron born in fission appears in energy group :math:`g` — carries a
**probability-simplex law** for a fissile material, and is identically
zero (the **null spectrum**) for a non-fissile material:

.. math::
   :label: emission-spectrum-simplex

   \text{producing }(\nu\Sigma_f>0):\quad \sum_g \chi_g = 1,\ \ \chi_g \ge 0\ \forall g;
   \qquad
   \text{non-producing:}\quad \chi_g \equiv 0\ \forall g.

.. vv-status: emission-spectrum-simplex documented
.. (vv-status rationale) Definitional probability identity: the simplex / null
   law of the fission emission spectrum (Σχ=1, χ≥0 for producing; χ≡0 for
   non-producing). A probability-distribution invariant carried by the
   EmissionSpectrum type, gated by the foundation tests
   ``tests/data/test_emission_spectrum.py`` (``assert_normalized`` positive /
   negative legs + ``assert_null``). A definitional law, not a solver claim.

ORPHEUS encodes this law in
:class:`~orpheus.data.emission_spectrum.EmissionSpectrum`, a
:class:`numpy.ndarray` subclass that carries the law as inspectable
methods (``assert_normalized`` / ``is_emitting`` / ``assert_null``) — so
``mixture.chi.assert_normalized()`` reads like a law the data carries.
Both :class:`~orpheus.data.micro_xs.isotope.Isotope` and
:class:`~orpheus.data.macro_xs.mixture.Mixture` coerce ``chi`` to this
type at construction (``__post_init__``) via the single-source helper
:func:`~orpheus.data.emission_spectrum.enforce_emission_spectrum`, which
cross-checks the law conditionally on ``is_producing`` — keyed on
PRODUCTION (:math:`\nu\Sigma_f > 0`), because :math:`\chi` is consumed
only as a fission *source* :math:`\chi\,\nu\Sigma_f\,\varphi`, so a valid
simplex is required exactly where production is nonzero. Enforced ONCE,
at the data source, NEVER per-cell (a non-fissile cell legitimately holds
:math:`\chi = 0` under per-cell broadcast). The value-object is
deliberately NOT a member of the flux/operator algebra: it is a validated
newtype, not a per-cell field type.

Why the law is enforced at the source (and never per-cell)
----------------------------------------------------------

The simplex / null law of :eq:`emission-spectrum-simplex` is a property
of a **material**, not of a **cell**. ORPHEUS therefore enforces it
exactly once — at the point where a material's cross sections become a
data object — and never again downstream.

There are exactly two places where a material is born:
:meth:`~orpheus.data.micro_xs.isotope.Isotope.__post_init__` (a single
isotope from a GENDF library) and
:meth:`~orpheus.data.macro_xs.mixture.Mixture.__post_init__` (a
homogenized mixture). Both route ``self.chi`` through the same
single-source law
:func:`~orpheus.data.emission_spectrum.enforce_emission_spectrum`, which
coerces the bare array to an
:class:`~orpheus.data.emission_spectrum.EmissionSpectrum` and runs the
conditional cross-check:

.. code-block:: python

   def enforce_emission_spectrum(chi, *, is_producing):
       chi = EmissionSpectrum(chi)
       chi.assert_normalized() if is_producing else chi.assert_null()
       return chi

The two ``__post_init__`` bodies were duplicated before S10a; collapsing
them into one helper (single source of truth) is what guarantees the
isotope guard and the mixture guard can never drift apart — they are
literally the same three lines.

**Why not per-cell?** Every transport solver eventually broadcasts the
material spectrum across the spatial mesh — ``chi[None, :]`` in the
fission-source assembly, where the leading axis is the cell index. In a
heterogeneous geometry, the cells that fall in a non-fissile region
legitimately carry :math:`\chi = 0`. A per-cell validator would have to
permit :math:`\chi = 0` *everywhere*, which makes the simplex clause
vacuous — it could no longer distinguish a genuine null spectrum (a
non-producing material) from a malformed one (a producing material whose
spectrum was accidentally zeroed). The discriminating information lives
upstream, at the material: only there do we know whether *this*
material is supposed to emit. Validating at the broadcast site discards
exactly the predicate (``is_producing``) that makes the law meaningful.

This is the conceptual boundary between
:class:`~orpheus.data.emission_spectrum.EmissionSpectrum` and the
(#263-deferred) per-cell ``SpectrumField``. ``EmissionSpectrum`` is a
**validated value-object**: a thin :class:`numpy.ndarray` subclass that
carries the material-level law as inspectable methods, with no membership
in the flux/operator algebra (no ``__add__`` gate, no cross-mesh partner
check, no change-of-basis morphism). It validates the values *and lets
them ride unchanged* through every existing call site, because an
``EmissionSpectrum`` *is* an ndarray — ``chi[None, :]``, ``chi.sum()``,
``chi.copy()``, and einsum-feeding the fission source all keep working
with zero ripple. A per-cell ``SpectrumField``, by contrast, *would*
live in the field-typed algebra, and its broadcast of :math:`\chi`
across non-fissile cells is precisely where the simplex law stops
applying. The law belongs to the material; the field would belong to
the mesh.


The tolerance rationale: atol = 1e-12 over a real FP residual
-------------------------------------------------------------

:meth:`~orpheus.data.emission_spectrum.EmissionSpectrum.assert_normalized`
checks the sum clause with
``np.isclose(self.sum(), 1.0, atol=1e-12, rtol=0)``. The absolute
tolerance is set by a measurement, not a guess.

Real GENDF data is not exactly normalized in IEEE-754. Probing the
shipped library directly:

.. code-block:: python

   >>> from orpheus.data.micro_xs import load_isotope
   >>> iso = load_isotope("U_235", 294)
   >>> iso.chi.sum()
   1.0000000000000002          # FP residual ≈ 2.22e-16 (1 ULP at 1.0)

The summed spectrum lands one ULP above unity — a residual of
:math:`\approx 2.22\times10^{-16}`. The tolerance must clear this FP
noise comfortably while still failing on a *physical* normalization
error (a transcription bug, a wrong group count, a renormalization that
forgot a group). The band is chosen as:

.. list-table:: Choosing the normalization tolerance
   :header-rows: 1
   :widths: 32 22 22

   * - Scale
     - Magnitude
     - Verdict
   * - Real GENDF FP residual (U-235)
     - :math:`\approx 2.2\times10^{-16}`
     - MUST pass
   * - ``atol = 1e-12`` (chosen)
     - :math:`10^{-12}`
     - the threshold
   * - A :math:`10^{-6}`-scale physical error
     - :math:`10^{-6}`
     - MUST fail

The chosen ``atol`` sits roughly four orders of magnitude above the
observed FP residual and six orders below a realistic physical error,
so it is robust at both ends. The ``rtol=0`` choice is deliberate: at
:math:`\sum_g \chi_g \approx 1`, an absolute tolerance is the natural
metric (the target value is exactly 1.0, so relative and absolute
coincide to leading order; pinning ``rtol=0`` removes any ambiguity).
Both ends of this band are pinned by foundation tests
(``test_one_ulp_residual_passes`` mirrors the 2.22e-16 residual and must
NOT raise; ``test_off_by_1e6_raises`` injects the :math:`10^{-6}` error
and must raise) in ``tests/data/test_emission_spectrum.py``.

The **null** spectrum is held to a stricter standard.
:meth:`~orpheus.data.emission_spectrum.EmissionSpectrum.assert_null`
checks ``np.all(self == 0)`` — *exact* zero, no tolerance. A null
spectrum is never the result of a summation that might accumulate FP
error; it is a *constructed* zero (``np.zeros(ng)``). There is no
physical residual to absorb, so demanding exact zero catches any
accidental nonzero entry (even a :math:`10^{-9}` leak) that a tolerant
check would silently pass.

The two clauses of ``assert_normalized`` — the sum clause and the
non-negativity clause — are checked **independently**. A spectrum like
``[1.2, -0.2]`` sums to exactly 1.0 yet is not a valid probability
distribution; checking only the sum would pass it. The independent
``np.all(self >= 0)`` clause catches it on the negativity test, with a
distinct error message naming the negative entries. (Foundation leg
``test_negative_entry_raises_even_when_sum_is_one``.)


.. _emission-spectrum-production-keying:

Why the law keys on production, not fissionability
--------------------------------------------------

The conditional in ``enforce_emission_spectrum`` branches on
``is_producing`` (:math:`\nu\Sigma_f > 0`), **not** on fissionability
(:math:`\Sigma_f > 0`). This is the load-bearing design decision of
S10a, and it follows directly from *how* :math:`\chi` is consumed.

The emission spectrum appears in the transport equation only inside the
fission **source** — the rate at which fission neutrons are born into
group :math:`g`:

.. math::
   :label: emission-spectrum-fission-source

   q^{\mathrm{fiss}}_g(\mathbf{r}) \;=\; \chi_g \sum_{g'} \nu\Sigma_{f,g'}\,\varphi_{g'}(\mathbf{r}).

.. vv-status: emission-spectrum-fission-source documented
.. (vv-status rationale) Governing / definitional equation: the fission source
   term (χ_g times the production rate) — the standard transport definition that
   makes the emission-spectrum law key on production (νΣf), not fissionability.
   Realised by FissionOperator and exercised in the multi-group eigenvalue
   chain; the production keying is gated by the S10a foundation suite. A
   governing definition, not a separate solver claim.

:math:`\chi` is never used on its own; it is always multiplied by the
**production** cross section :math:`\nu\Sigma_f`. Therefore the only
quantity that determines whether :math:`\chi` matters is precisely
:math:`\nu\Sigma_f`: where :math:`\nu\Sigma_f > 0`, the spectrum must be
a genuine probability distribution (the birth density must integrate to
the production rate); where :math:`\nu\Sigma_f = 0`, the product
:math:`\chi\,\nu\Sigma_f` vanishes regardless of :math:`\chi`, so the
law that has teeth is precisely the null spectrum. **A valid emission
spectrum is required exactly where production is nonzero** — which is
the statement of :eq:`emission-spectrum-simplex` with the
``is_producing`` predicate as its keying.

Fissionability is a *different* question. A material is fissile
(:math:`\Sigma_f > 0`) if it *can* undergo fission; it is producing
(:math:`\nu\Sigma_f > 0`) if it *emits* neutrons when it does. ORPHEUS
keeps both predicates, with distinct roles:

.. list-table:: The two fission predicates
   :header-rows: 1
   :widths: 22 26 30

   * - Predicate
     - Definition
     - Role
   * - ``is_fissile``
     - :math:`\Sigma_f > 0` (``np.any(sigF > 0)``)
     - "Can it fission?" — consumed by
       :func:`~orpheus.data.macro_xs.mixture.compute_macro_xs` to
       auto-detect ``fissile_indices`` (which isotopes contribute to the
       homogenized production XS).
   * - ``is_producing``
     - :math:`\nu\Sigma_f > 0` (``np.any(nubar * sigF > 0)`` /
       ``np.any(SigP > 0)``)
     - "Does it emit fission neutrons?" — the predicate the
       :math:`\chi` simplex/null law keys on.

For **real** nuclear data flowing through
:func:`~orpheus.data.macro_xs.mixture.compute_macro_xs`, the two
predicates coincide: the production XS is computed as
:math:`\Sigma_P = \overline{\nu}\,\Sigma_f` with
:math:`\overline{\nu} > 0` for every fissile nuclide, so
:math:`\nu\Sigma_f > 0 \;\Leftrightarrow\; \Sigma_f > 0`. On the real
GENDF path, keying on either predicate gives the same answer. This is
proven on shipped data by ``TestRealGendfConstructs`` in
``tests/data/test_chi_invariant_enforcement.py``: U-235 reports
``is_producing`` and its 2.22e-16-residual spectrum passes the simplex
clause; O-016 and H-001 report *not* ``is_producing`` and their
genuinely-zero spectra pass the null clause.

The two predicates **diverge** for a synthetic fixture that sets the
production XS directly while leaving the fission XS at zero — a
*production-bearing but non-fissile* material. The canonical instance is
the trajectory-resolvent billiard reference solver
(``tests/derivations/test_trajectory_resolvent_billiard.py``), which
builds a multiplying medium by setting ``SigP = νΣ_f > 0`` and
``SigF = 0`` (the billiard solver reads ``SigP`` and ``chi`` to build
the fission source; it never reads ``SigF``). For this fixture:

- ``is_fissile`` (keyed on ``SigF``) would return **False**, demanding a
  null spectrum — yet :math:`\chi` is genuinely consumed via
  :math:`\chi\cdot\mathtt{SigP}`, so zeroing it would corrupt the
  multi-group billiard result. Keying the guard on fissionability would
  force a *behavior-changing* edit to a verification reference.
- ``is_producing`` (keyed on ``SigP``) returns **True**, correctly
  requiring the simplex law. The fixture's default simplex passes the
  guard with no special handling.

Production-keying classifies the production-bearing fixture correctly
**with no stand-in**. An earlier draft of S10a keyed the guard on
``is_fissile`` and worked around this divergence by patching the billiard
fixture to carry a synthetic ``SigF`` (a stand-in for the untracked true
:math:`\Sigma_f = \mathtt{SigP}/\nu`). Re-keying on production *retired*
that hack: the fixture now sets ``SigF = 0`` honestly, and the
divergence that looked like an irreducible design seam dissolves. The
isolation is pinned directly by ``test_non_producing_nonzero_raises``,
which builds an isotope with ``sigF > 0`` but ``nubar = 0`` (fissile, yet
non-producing) and asserts its spectrum must be null — a leg that
*requires* ``is_producing`` to differ from ``is_fissile``.

The lesson generalizes beyond :math:`\chi`: a guard-at-the-source should
key on the predicate that matches **why the guarded value exists**. The
emission spectrum exists to be a fission source, so production is the
right key. Fissionability — which the macro-XS pipeline still needs for
``fissile_indices`` — is the wrong key for the :math:`\chi` law, and
choosing it manufactures a seam where there is none.


The behavior-neutral precursor: zeroing dead non-fissile χ
----------------------------------------------------------

Turning on a strict source-level guard is only safe if every existing
call site already satisfies it. An audit before S10a found that the
:ref:`synthetic verification cross-section library <synthetic-xs-library>`
(regions A/B/C/D) carried **dead** nonzero spectra: the non-fissile
regions B, C, and D held placeholder spectra such as ``chi=[1.0]``,
``chi=[1, 0]``, and ``chi=[0.6, 0.35, 0.05, 0]`` — left over from an
era when every region was given a spectrum by default. With the
production-keyed guard active, these non-producing regions
(:math:`\nu\Sigma_f = 0`) would demand the *null* spectrum and the guard
would reject them on construction, reddening the SN / CP / MoC suites en
masse.

The precursor zeroed every such dead spectrum (region A, the only
producing region, is untouched). This change is **behavior-neutral by
construction**: a non-producing material has :math:`\nu\Sigma_f \equiv 0`,
so its contribution to the fission source
:eq:`emission-spectrum-fission-source` is
:math:`\chi\,\nu\Sigma_f\,\varphi \equiv 0` *no matter what* :math:`\chi`
holds. The spectrum of a non-producing region is multiplied by zero
everywhere it is consumed; zeroing :math:`\chi` itself merely makes the
dead value explicit. Nothing observable depends on it.

This inertness is not asserted — it is **proven** by the byte-identical
diamond-difference regression in
``tests/sn/regression/test_dd_regression.py``. That gate pins the SN
flux/eigenvalue snapshots with an exact-equality contract. After zeroing
the non-fissile spectra, every snapshot reproduces bit-for-bit (no
snapshot moved). A moved snapshot would have meant :math:`\chi` was *not*
inert — that some non-producing region's spectrum was leaking into a
result — and the precursor would have been wrong. None moved, so the
zeroing is empirically confirmed behavior-neutral, exactly the standard
of evidence the bit-identity-vs-principled-equivalence discipline
demands: a regression contract is the right gate when the change is
supposed to alter *nothing*.

The audit also taught a scope lesson worth recording: a
guard-at-the-data-source has a blast radius equal to **every direct
constructor** of the guarded type, not just the factory path. Regions
reached through ``get_mixture("B"/"C"/"D")`` were fixed transitively by
zeroing the library, but a further set of test fixtures called
``Mixture(...)`` / ``make_mixture(...)`` *directly* with inline
placeholder spectra — invisible to a factory-path audit. A deterministic
static scan of every constructor call (classifying the production-keying
field against the guarded :math:`\chi` field) found them all in one pass,
where running the suite would have leaked them one failure at a time.


The S10b seam: production-weighted multi-fissile χ_mix
------------------------------------------------------

S10a installs the type and the source-level guard but changes **no**
:math:`\chi` value. In particular,
:func:`~orpheus.data.macro_xs.mixture.compute_macro_xs` still homogenizes
a multi-isotope mixture by taking the spectrum of the **first** fissile
isotope verbatim:

.. code-block:: python

   # Fission spectrum — use first fissile isotope's chi (simplification)
   chi = np.zeros(NG)
   if fissile_indices:
       chi = isotopes[fissile_indices[0]].chi.copy()

This is correct for a single-fissile mixture but a known simplification
when two or more fissile isotopes coexist (e.g. a U-235 / Pu-239 fuel),
since each has a distinct birth spectrum. The physically correct mixture
spectrum is the **production-weighted average** of the per-isotope
spectra:

.. math::
   :label: emission-spectrum-chi-mix

   \chi^{\mathrm{mix}}_g
   \;=\;
   \frac{\displaystyle\sum_{i\in\text{producing}} \chi^{(i)}_g \sum_{g'} N_i\,\nu\sigma^{(i)}_{f,g'}\,\varphi_{g'}}
        {\displaystyle\sum_{i\in\text{producing}} \sum_{g'} N_i\,\nu\sigma^{(i)}_{f,g'}\,\varphi_{g'}},

a convex combination of the per-isotope spectra with weights equal to
each isotope's share of the total production rate. (In practice the
spectrum is condensed against a representative flux or assumed
group-independent in weight; the structure is a convex combination
regardless.)

S10a is the deliberate seam at which S10b will attach. The decisive
property is that **a convex combination of probability simplices is
itself a probability simplex**: if each :math:`\chi^{(i)}` satisfies
:math:`\sum_g \chi^{(i)}_g = 1` with :math:`\chi^{(i)}_g \ge 0`, and the
weights :math:`w_i \ge 0` sum to 1, then
:math:`\sum_g \sum_i w_i \chi^{(i)}_g = \sum_i w_i = 1` and every entry
is a non-negative sum of non-negative terms. Therefore
:math:`\chi^{\mathrm{mix}}` will pass
:meth:`~orpheus.data.emission_spectrum.EmissionSpectrum.assert_normalized`
**verbatim** — the simplex validator written for S10a is already the
correctness floor for the S10b weighting, with no new law required. The
value-object was designed so that the harder multi-fissile computation
inherits its invariant for free; S10b implements the weighting and
trusts the same guard to catch a mis-weighted result (one whose weights
do not sum to 1, or that admits a negative entry).


.. _emission-spectrum-verification:

Verification chain
------------------

The simplex / null law is a software invariant — a property of a
probability distribution, not a theory-page equation about a solver —
so the gates are ``@pytest.mark.foundation`` and carry no
``verifies(...)`` edge (there is no flux or eigenvalue claim to pin).

.. list-table:: Foundation gates for the emission-spectrum law
   :header-rows: 1
   :widths: 38 40

   * - Gate
     - What it pins
   * - ``tests/data/test_emission_spectrum.py``
     - The **intrinsic** simplex / null law on the value-object itself
       (vv anti-pattern #11, both legs): each validator gets a positive
       leg (correct instance MUST NOT raise) and a negative leg (broken
       instance MUST raise), with simplex / non-simplex arrays built by
       hand (L11 structural independence — never via the production
       builder the guard protects). Includes the tolerance-band legs
       and the ndarray-subclass zero-ripple legs.
   * - ``tests/data/test_chi_invariant_enforcement.py``
     - The **container** cross-check: for each of ``Mixture`` and
       ``Isotope``, both branches (producing → simplex, non-producing →
       null) and both legs, plus ``test_non_producing_nonzero_raises``
       (the ``is_producing`` ≠ ``is_fissile`` isolation leg) and the
       real-GENDF ``TestRealGendfConstructs`` no-red proof on U-235 /
       O-016 / H-001.

Both gates run under the canonical ``python -O`` invocation. Every
assertion routes through ``pytest.raises`` / ``np.testing.*`` (or the
``_require`` helper that wraps ``pytest.fail``) rather than a bare
``assert`` — so no leg is silently stripped under ``-O`` (vv failure
mode 8, the compiled-out assertion). The byte-identical DD regression
(``tests/sn/regression/test_dd_regression.py``) is the third leg: it
proves the precursor zeroing of dead non-fissile :math:`\chi` changed no
solver output.


.. seealso::

   - :ref:`theory-homogeneous` — first consumer of the XS pipeline;
     demonstrates the full path from ``load_isotope()`` to :math:`k_\infty`.
   - :ref:`theory-verification` — verification uses :ref:`synthetic cross
     sections <synthetic-xs-library>` (regions A/B/C/D), not this pipeline.
   - :ref:`theory-collision-probability`, :ref:`theory-discrete-ordinates`,
     :ref:`theory-method-of-characteristics`, :ref:`theory-monte-carlo` — all
     transport solvers consume ``Mixture`` objects from this pipeline.
