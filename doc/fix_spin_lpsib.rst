.. index:: fix spin/lpsib

fix spin/lpsib command
======================

Syntax
------

.. parsed-literal::

   fix ID group-ID spin/lpsib T alpha seed lattice frozen|moving accept N [options]

Options are ``predictor lagged|fresh``, ``ou T damp seed``,
``box iso P W`` and ``box_ou T damp seed``. Seeds for the spin, atomic OU and
box OU streams must be distinct. Temperatures are in K, damping times in ps,
pressure P in bar and W is the inertia of the volume coordinate in metal
units (eV ps^2/Angstrom^6). W is positive; it is not a relaxation time.
The default ``predictor lagged`` uses the cached field. ``predictor fresh`` selects centered SIB.
Frozen fresh is a splitting control and timing mode.

Description
-----------

This fix implements the cached-field integrator (CFI). CFI reuses the
previous normalized midpoint field for prediction and
makes one core pair evaluation per step. Centered SIB obtains a fresh
predictor field and makes two. The first run step uses the setup field.
``accept N`` requests an additional accepted-state evaluation every N steps;
zero disables it. This output evaluation does not enter the dynamics.
To read energy or conservative spin forces, align dump and accept intervals.

Frozen mode changes spins only. Moving mode integrates positions and velocities
with centered spin-lattice forces. Use the built-in atomic OU half-steps for
lattice thermalization; an external lattice thermostat is refused. Conservative
``precession/spin`` fields are included, including with hybrid spin pairs.
The box option uses an isotropic canonical-volume Andersen coordinate; box OU
adds an independent OU bath for its conjugate momentum. The implementation
keeps the box shape fixed and requires an orthogonal periodic 3D cell.

Restrictions
------------

Exactly one MPI rank, Verlet, metal units and atomic spins are required.
Moving mode requires group all and newton on. Do not overlap time integrators,
spin thermostats or box-changing fixes. The NPT neighbor policy must be
``neigh_modify every 1 delay 0 check yes``. NPT admits only ``zero``, ``lj/cut``,
``spin/exchange`` and their hybrid combinations. ML NPT is rejected by the
pair-virial guard. NequIP's newton-off requirement is incompatible with moving
CFI. NPT has been tested only for analytic Fe, not ML potentials.

Deterministic second order refers to the undamped, unthermostatted equations.
No stochastic weak-order claim is made. Alpha equal to one is used for
thermodynamic sampling only. See :doc:`SPINGNN` and the repository limitations.

Output
------

The first three fix vector entries are completed steps, accepted-state
output evaluations and maximum spin-norm deviation. Additional box diagnostics
are defined in ``FixSpinLPSIB::compute_vector``. A run prints a separate
``OBSERVE85`` call ledger; see the repository VALIDATION guide.
