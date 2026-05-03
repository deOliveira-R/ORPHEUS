# `bn_method/` — B_N method — reserved (lower priority)

**Status:** empty / placeholder. Implementation deferred. Lower
priority than the other reserved folders.

## Intent

The B_N method is the sister of F_N: it collocates the F_N equations
on the **boundary points** instead of interior collocation points.
Used historically for diffuse-reflection benchmarks and for problems
where the boundary structure dominates the angular response.

In the ORPHEUS verification stack, F_N already covers most of the
slab/sphere benchmark surface — so B_N is a niche addition rather
than a foundational pillar.

## Canonical references (LOCAL in `scratch/literature/`)

* **Brockmann, R. (1981).** B_N method for anisotropic-scattering
  treatment. Local PDF: `scratch/literature/Brockmann1981.pdf`.

Additional context: Sood/Forster/Parsons 2003 cite the B_N method
in the LA-13511 benchmark catalogue context. The Garibba–Rojas
1980s technical reports (NOT local) contain the original
diffuse-reflection benchmark formulations.
