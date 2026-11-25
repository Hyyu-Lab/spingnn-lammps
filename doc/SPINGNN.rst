SPINGNN package
===============

This overlay contains the cached-field integrator (CFI, fix style
``spin/lpsib``) with frozen, moving-lattice and isotropic
canonical-volume modes, the centered SIB control and spin-aware ML pair
relocation. Native ``spin/sib`` and ``spin/sib_npt`` remain in SPIN.
Both ``PKG_SPIN`` and ``PKG_SPINGNN`` must be enabled in CMake.

.. toctree::
   :maxdepth: 1

   fix_spin_lpsib
   fix_spin_sib
   fix_spin_sib_npt
   pair_spin_allegro_spin
   pair_spin_nequip_spin

Use the examples as demonstrations, not timestep recommendations. The offline
geometric-probe spin temperature (GPST), T_GP, lives in ``tools/gpst`` in the
overlay repository. Its accepted-state dump adapter requires an explicit hbar
conversion and checks unit spin directions. Model files are not distributed.
See README, VALIDATION and KNOWN-LIMITATIONS for the tested support boundary.
