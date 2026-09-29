r"""F_N method analytical benchmark family — Sood/Forster/Parsons (2003).

This package is a structurally-independent verification ground for
ORPHEUS transport solvers, mounted on the published Sood/Forster/Parsons
analytical benchmark test set :cite:`SoodForsterParsons2003`. The benchmark
catalogue contains 75 critical configurations (43 1G, 30 2G, 1 3G, 1 6G;
24 infinite, 24 slabs, 9 cylinders, 14 spheres, 4 ISLC) whose critical
radii, :math:`k_\infty`, and selected scalar-flux ratios were computed
in the peer-reviewed transport-theory literature using **Case singular
eigenfunctions, F_N method, B_N method, and Green's-function method**.

The package follows the algebra-of-record bifurcation discipline:

* :mod:`.origins` — Branch 1 SymPy: closed-form derivations from the
  transport equation. These verify the *algebra* — that the published
  closed forms follow from the reduction chain documented in Sood 2003
  Section 3 + Appendix A.
* :mod:`.multi_group` — Branch 2 numpy/scipy: production code that
  evaluates the closed forms on the catalogued cases. These verify the
  *numerics* — that the algebra evaluates to the published reference
  values.
* The machine-readable catalogue of the Sood cases is
  :mod:`orpheus.derivations.continuous.sood_registry`.

The F_N solvers, built from the primary papers Sood cites (the Sood
paper states the truth set, not the method):

* :mod:`.core` — the shared F_N primitives: the Case dispersion-relation
  roots, the F_N collocation matrix and the half-range moment integrals.
* :mod:`.slab` — the one-group bare critical slab (Siewert-Benoist 1979,
  Grandjean-Siewert 1979), its interior flux (the Kaper-Lindeman-Leaf
  1974 recipe), and the one-group reflected slab (Neshat-Maiorino 1980).
* :mod:`.sphere` — the one-group bare critical sphere (Siewert-Thomas
  1986) and its interior flux.
* :mod:`.cylinder` — a placeholder. The bare cylinders' critical
  dimensions come from
  :mod:`orpheus.derivations.continuous.singular_eigenfunction.cylinder`
  (Westfall-Metcalf 1972).
* :mod:`.moment_space` — :class:`~.moment_space.MomentSpace`, the
  reference generator that reads a ``StructuredGeometry`` and routes it
  to these solvers.

"""
from __future__ import annotations
