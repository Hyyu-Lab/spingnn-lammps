.. index:: fix spin/sib

fix spin/sib command
====================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID spin/sib T alpha seed lattice moving|frozen

Example
"""""""

.. code-block:: LAMMPS

   fix integrate all spin/sib 50 1 531791 lattice moving

Description and restrictions
""""""""""""""""""""""""""""

T is the spin target in K; alpha is dimensionless Gilbert damping. The seed selects the spin RNG. Zero T and alpha provide the unthermostatted spin path. Moving integrates atomic positions/velocities; frozen does not. Do not overlap lattice integration fixes. This is the normalized predictor-midpoint variant, not Mentink Eq. (18); order theorems do not transfer automatically. The inherited implementation lives in SPIN/fix_spin_sib.cpp. See the analytic NVE and lattice-Langevin NVT examples.

Source anchor: ``SPIN/fix_spin_sib.cpp:80-143``.

All scientific/support boundaries are mandatory: :doc:`SPINGNN` links the
package limitations, examples, build and provenance notes.
