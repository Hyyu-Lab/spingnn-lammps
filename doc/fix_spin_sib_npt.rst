.. index:: fix spin/sib_npt

fix spin/sib_npt command
========================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID spin/sib_npt temp Tstart Tstop Tdamp pressure-keywords alpha A seed S

Example
"""""""

.. code-block:: LAMMPS

   fix integrate all spin/sib_npt temp 50 50 0.1 iso 0 0 1 alpha 1 seed 531791

Description and restrictions
""""""""""""""""""""""""""""

Uses the normal FixNH pressure forms and integrates lattice plus stochastic spins. In metal units Tdamp/Pdamp are ps and pressure is bar. The spin target is Tstop, not an independent bath schedule. Alpha-first and seed-first both parse. Do not apply another time integrator to the same atoms. The original cached-field behavior is preserved. The former opt-in cache-refresh style from earlier development is not included in this release. Centered SIB is selected only through ``spin/lpsib predictor fresh``; it does not change this native fix.

Source anchor: ``SPIN/fix_spin_sib_npt.cpp:80-130``.

All scientific/support boundaries are mandatory: :doc:`SPINGNN` links the
package limitations, examples, build and provenance notes.
