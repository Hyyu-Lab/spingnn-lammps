.. index:: pair_style spin/allegro_spin

pair_style spin/allegro_spin command
====================================

Syntax
""""""

.. parsed-literal::

   pair_style spin/allegro_spin
   pair_coeff * * model-file element ...

Example
"""""""

.. code-block:: LAMMPS

   pair_style spin/allegro_spin
   pair_coeff * * external-model.pt2 Cr S Br

Description and restrictions
""""""""""""""""""""""""""""

Requires spin atoms, newton pair on and compatible deployed metadata. Element names must match atom types. Model paths are user supplied; no trained model is distributed here. TorchScript/PT2 support and runtime-device validation remain model/deployment specific. The Kokkos variant spin/allegro_spin/kk remains in KOKKOS, conditioned on SPINGNN. fm uses the unit-direction angular-frequency convention in the package limitations. Do not use Nurdin as an ML temperature estimator. No GPU/ML regression claim follows from CPU registration or parse checks.

Source anchor: ``SPINGNN/pair_allegro_spin.cpp:145-198``.

All scientific/support boundaries are mandatory: :doc:`SPINGNN` links the
package limitations, examples, build and provenance notes.
