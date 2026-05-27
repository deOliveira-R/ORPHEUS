Architecture
============

ORPHEUS organizes its code by **mathematical knowledge layer** rather than
by feature or by file size. A module's home is the lowest-knowledge layer
whose vocabulary suffices to define it; imports flow only from
more-knowledge to less-knowledge. This chapter records the criterion, the
package-to-layer assignment, and the import-linter test that enforces it.

.. toctree::
   :maxdepth: 2

   layering
