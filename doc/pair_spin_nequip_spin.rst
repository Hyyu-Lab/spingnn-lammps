.. index:: pair_style spin/nequip_spin

pair_style spin/nequip_spin command
===================================

Syntax
""""""

.. parsed-literal::

   pair_style spin/nequip_spin
   pair_coeff * * model-file element ...

Example
"""""""

.. code-block:: LAMMPS

   pair_style spin/nequip_spin
   pair_coeff * * external-model.pth Fe O

Description and restrictions
""""""""""""""""""""""""""""

Requires spin atoms, newton pair off and compatible deployed NequIP metadata. Element names map every atom type. Supply the deployed model externally; there are no model weights in these examples. Metadata, neighborhood and output validation remain in the original implementation. Use the package force conventions; no general Nurdin/ML thermometer is valid. CPU smoke checks do not certify GPU/ML equivalence.

Source anchor: ``SPINGNN/pair_nequip_spin.cpp:105-160``.

All scientific/support boundaries are mandatory: :doc:`SPINGNN` links the
package limitations, examples, build and provenance notes.
