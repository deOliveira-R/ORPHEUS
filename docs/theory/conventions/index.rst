.. _theory-conventions:

===========
Conventions
===========

**Read this part first.** It is the one part with no counterpart in the
literature, and the one whose absence produces the worst failure class in a
transport code: *code that runs, converges, and is wrong by a constant.*

A conventions crosswalk is not pedantry here — the published canon does not
agree with itself, and largely does not warn you:

- **Quadrature weights.** Hébert §3.9.1 normalizes 1-D Gauss-Legendre to
  :math:`\sum w = 2`; five pages later, Eqs. (3.363)–(3.364) normalize to
  :math:`\sum w = 1` over the positive octant — same symbol :math:`w_n`, two
  normalizations, no note. ORPHEUS's ERR-025 was exactly a missing
  :math:`1/W`, masked because Gauss-Legendre's :math:`W = 2` made the wrong
  formula accidentally right.
- **The** :math:`(2\ell+1)` **prefactor** carries :math:`4\pi` in Hébert
  (3.30) but :math:`2` in (3.425) — the same object, silently tied to
  dimensionality. ORPHEUS's ERR-039 was this prefactor.
- **The scattering arrow points opposite ways across texts.** Hébert writes
  :math:`\Sigma_s(E \leftarrow E', \Omega \leftarrow \Omega')`
  (destination-first); Bell & Glasstone write :math:`\sigma(x, E' \to E)` and
  Stacey :math:`\Sigma_s(\Omega' \to \Omega)` (source-first).

Each is a live convention trap, and each has a catcher in the ERR catalog.
Every import of an external equation must cross the boundary through this
part.

.. note::

   Still growing: the normalization page (weight sums, the
   :math:`(2\ell+1)` prefactor family, the :math:`\alpha`-recursion
   crosswalk) and the split of the indexing / cross-section page are
   queued — see issue
   `#231 <https://github.com/deOliveira-R/ORPHEUS/issues/231>`_.

.. toctree::
   :maxdepth: 2

   notation
   index_convention
