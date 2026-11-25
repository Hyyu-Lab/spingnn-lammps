# SpinGNN-LAMMPS

[![arXiv](https://img.shields.io/badge/arXiv-2611.32129-b31b1b.svg)](https://arxiv.org/abs/2611.32129)
[![License: GPL-2.0-only](https://img.shields.io/badge/license-GPL--2.0--only-blue.svg)](LICENSE)
[![GPST tools: AGPL-3.0-only](https://img.shields.io/badge/GPST_tools-AGPL--3.0--only-blue.svg)](tools/gpst/LICENSE)

A LAMMPS package for finite-temperature spin–lattice dynamics with
machine-learned and analytic magnetic potentials. It is the code for

> Y. Zheng, L. Hong, X. Gong, H. Xiang and H. Yu, *Efficient spin–lattice
> dynamics and spin thermometry with machine-learned magnetic potentials*,
> [arXiv:2611.32129](https://arxiv.org/abs/2611.32129) (2026).

The package provides two methods:

| Method | What it does | Model evaluations |
|---|---|---|
| **CFI**, the cached-field integrator (LAMMPS `fix spin/lpsib`) | Integrates coupled spin and lattice dynamics on a frozen lattice, on a moving lattice (NVE, NVT) and at constant pressure (isotropic NPT). The spin predictor reuses the field cached from the previous step. Second-order accurate for deterministic NVE. | One core evaluation per step; native SIB needs two |
| **GPST**, the geometric-probe spin temperature (`tools/gpst/`) | Measures the spin temperature offline from stored unit spins and conservative spin forces. The underlying identity holds exactly for any smooth Hamiltonian, including energies that are nonlinear in each spin. | None |

Both methods use only first derivatives of the energy: atomic forces, the
spin force and, for NPT, the virial. Neither uses a Hessian, a
Hessian–vector product or finite differences, and neither contains element-,
material- or model-specific code.

This repository is an add-on for [SpinGNN](https://github.com/Hyyu-Lab/SpinGNN).
An installer copies the package into a SpinGNN checkout and applies a small
patch to its LAMMPS. No model files are distributed.

## Contents

| Path | Contents |
|---|---|
| `src/SPINGNN/` | Package source: `fix spin/lpsib` (CFI), package notes and known limitations |
| `patches/lammps.patch` | Changes to LAMMPS core files: package registration for CMake and make, model-call counters in Verlet and native SIB, the lattice-force switch used by moving-lattice CFI, a spin-mapping fix in `set`, and a modification notice |
| `cmake/SPINGNN.cmake` | CMake module for the package |
| `doc/` | Manual pages for `fix spin/lpsib`, `fix spin/sib`, `fix spin/sib_npt` and the spin-aware Allegro and NequIP pair styles |
| `examples/SPINGNN/` | Input scripts for bcc Fe with analytic exchange (CFI and native SIB) and a CrSBr template for an Allegro model |
| `tools/gpst/` | GPST: reference and accelerated Python implementations, dump adapter, inference script and tests |
| `docs/USING-OTHER-POTENTIALS.md` | How to use your own potential |
| `install.sh` | Installer |
| `VALIDATION.md` | What has been tested, and how |

## Requirements

| Item | Requirement |
|---|---|
| SpinGNN | Tag [`v0.1.0`](https://github.com/Hyyu-Lab/SpinGNN/tree/v0.1.0) of Hyyu-Lab/SpinGNN, which provides LAMMPS with the SPIN package and the spin-aware ML pair styles. SpinGNN 2.x has a different layout and is not supported by this release. |
| Build | CMake, a C++17 compiler, MPI and a LibTorch build compatible with your models (SpinGNN depends on LibTorch) |
| Installer | Bash, Git and GNU coreutils (`sha256sum`) |
| GPST | Python 3 with NumPy and SciPy; pytest for the tests; optionally a C compiler for the native kernel |

The validated build used GCC 10.5.0 and OpenMPI 4.1.6 on Linux; see
[VALIDATION.md](VALIDATION.md).

## Installation

```sh
# 1. Get SpinGNN v0.1.0 and this package side by side
git clone --branch v0.1.0 https://github.com/Hyyu-Lab/SpinGNN.git SpinGNN
git clone https://github.com/Hyyu-Lab/spingnn-lammps.git

# 2. Install the package into the SpinGNN checkout
cd spingnn-lammps
./install.sh ../SpinGNN

# 3. Build LAMMPS with the SPIN and SPINGNN packages
cmake -S ../SpinGNN/lammps/cmake -B ../build \
  -D PKG_SPIN=ON -D PKG_SPINGNN=ON \
  -D CMAKE_PREFIX_PATH=/path/to/libtorch
cmake --build ../build --parallel
```

Before changing any file, `install.sh` checks the SHA-256 of every SpinGNN
file it modifies and of its own payload, refuses to overwrite existing files
and checks that the patch applies. It then moves the unchanged Allegro and
NequIP spin pair sources from SPIN into SPINGNN, copies the package files and
applies the patch. It refuses to install twice and has no rollback, so install
into a fresh checkout. The traditional make build is wired up but has not been
validated; use CMake.

## Quick start

The bcc Fe examples use an analytic potential and need no model:

```sh
cd ../SpinGNN/lammps/examples/SPINGNN
../../../../build/lmp -in in.fe-lpsib-moving-nvt
```

CFI runs on a single MPI rank. The fix line in that example is

```
fix integrate all spin/lpsib 300 1 531791 lattice moving accept 10 ou 300 .1 531792
```

| Part | Meaning |
|---|---|
| `300 1 531791` | Spin bath temperature (K), Gilbert damping α, seed of the spin noise |
| `lattice moving` | Integrate spins and atoms; `lattice frozen` integrates spins only |
| `accept 10` | Every 10 steps, evaluate the model once more at the accepted state for output. This evaluation does not enter the dynamics. Set it equal to the dump interval; `0` disables it. |
| `ou 300 .1 531792` | Built-in Ornstein–Uhlenbeck lattice bath: temperature (K), damping time (ps), seed |

Further options are `box iso P W` with `box_ou T damp seed` for isotropic NPT
(P in bar, W the inertia of the volume coordinate in eV ps²/Å⁶) and
`predictor fresh` for the centered-SIB control. All seeds must be distinct.
The full syntax is in [doc/fix_spin_lpsib.rst](doc/fix_spin_lpsib.rst).

| Example | What it runs |
|---|---|
| `in.fe-lpsib-frozen` | CFI on a frozen lattice |
| `in.fe-lpsib-moving-nve` | CFI on a moving lattice, NVE |
| `in.fe-lpsib-moving-nvt` | CFI on a moving lattice with the built-in lattice bath |
| `in.fe-lpsib-npt` | CFI at constant pressure (isotropic NPT) |
| `in.fe-lpsib-centered` | Centered-SIB control (`predictor fresh`) |
| `in.fe-sib-nve`, `in.fe-sib-nvt`, `in.fe-sib-npt` | Native SIB, for comparison |
| `in.crsbr-lpsib` | Template for an Allegro model of CrSBr; set `DATA_PATH` and `MODEL_PATH` |

The timesteps in the examples are demonstrations, not recommendations.

## Choosing an integrator

| Selection | Meaning |
|---|---|
| `spin/sib` | Native SIB, with its original numerical step order |
| `spin/sib_npt` | Native SIB at constant pressure, with its original cached-field behavior |
| `spin/lpsib ... lattice frozen` | CFI on a frozen lattice |
| `spin/lpsib ... lattice moving` | CFI on a moving lattice |
| `predictor fresh` | Centered SIB, a control that takes a fresh predictor field (also in frozen mode) |
| `box iso P W` | Isotropic constant-pressure dynamics, P in bar |
| `ou T damp seed`, `box_ou T damp seed` | Independent built-in baths for the atoms and the box |
| `spin/allegro_spin`, `spin/nequip_spin` | Spin-aware ML pair styles from SpinGNN, unchanged |

CFI makes one core model evaluation per step; centered SIB and native SIB make
two. With ML potentials this about halves the time per step. The fix style is
named `spin/lpsib` after the method's earlier working name, lagged-predictor
SIB.

## Measuring the spin temperature with GPST

The CFI examples write dumps with atom ID, position, unit spin, moment
magnitude and conservative spin force `fm`, in the column order that the GPST
adapter expects. Neither GPST script evaluates a model.

```sh
python tools/gpst/dump_to_moments.py trajectory.dump chain.npz \
  --hbar "$HBAR_EV_PS" --cutoff "$PROBE_CUTOFF_ANGSTROM" \
  --accept-every 10 --force-contract-confirmed
python tools/gpst/infer.py --kb "$KB_EV_K" --dt "$SAVED_INTERVAL_PS" \
  --seed-ids seed-a seed-b seed-c seed-d seed-e seed-f \
  --records a.npz b.npz c.npz d.npz e.npz f.npz --output estimate.json
```

The first command reduces one trajectory (one independent chain) to the
moments GPST needs. The second combines several independent chains into a
temperature estimate with its uncertainty. Pass ħ (LAMMPS `hplanck/(2π)`, in
eV ps) and k_B (in eV/K) explicitly for your unit system. `--dt` is the
interval between saved frames in ps, not the MD timestep. Unit spins and
`fm` must come from the same accepted state, so keep `accept` equal to the
dump interval. Report a temperature only when the returned status is
`IDENTIFIED`. Details, tests and the optional native kernel are described in
[tools/gpst/README.md](tools/gpst/README.md).

## Works with any magnetic potential

CFI reaches the potential only through the LAMMPS `PairSpin` interface, so
any `spin/*` pair style can drive it, for example `spin/exchange` or
`spin/allegro_spin`. GPST needs only unit spin directions and the
conservative spin force at the same stored state. The paper benchmarks
analytic bcc Fe potentials and a machine-learned CrSBr potential. To use
another potential, start from
[docs/USING-OTHER-POTENTIALS.md](docs/USING-OTHER-POTENTIALS.md) and check the
limitations below.

## Counting model evaluations

With screen output on, each Verlet run prints one `OBSERVE85` JSON line with
counts of pair, communication and neighbor operations. CFI also prints `S85_NATIVE_CONSTANTS`, and box mode
prints `S85N_UNITS`. These tag names are kept unchanged from the version used
in the paper. The counters wrap existing calls and add no model evaluations.
Read `setup`, `core`, `fresh`, `sib_mid` and `accept` separately and check
both `entered` and `completed`; setup and accepted-state output are not
per-step work. See [VALIDATION.md](VALIDATION.md).

## Limitations

1. CFI requires one MPI rank, metal units and the Verlet run style.
2. Moving-lattice CFI requires `group all` and `newton on`. Use its built-in
   lattice bath; an external `fix langevin` on the same atoms is refused.
   `spin/nequip_spin` requires `newton off` and therefore cannot be used with
   moving-lattice CFI.
3. NPT supports isotropic, orthogonal, periodic boxes. It has been tested only
   with analytic Fe potentials, and ML potentials are refused in NPT.
4. Second-order accuracy refers to deterministic NVE without damping or
   baths. No weak-order claim is made for the stochastic dynamics. α = 1 is
   meant for thermodynamic sampling.
5. The conservative spin force is f = ħ·fm for unit spin directions only.
6. The simplified Nurdin temperature of LAMMPS `compute spin` assumes an
   energy linear in each spin and is not valid for general ML potentials; use
   GPST instead.
7. GPST has a known platform-dependent result on one degenerate test input
   whose combined estimate is at floating-point roundoff scale. Do not use an
   `IDENTIFIED` result as a temperature when the point estimate is at
   roundoff scale, its standard error is zero, or the denominator interval is
   unresolved. See [tools/gpst/KNOWN-LIMITATIONS.md](tools/gpst/KNOWN-LIMITATIONS.md).
8. This packaging was validated on CPU. GPU execution, multiple MPI ranks,
   anisotropic boxes and the traditional make build are not certified here.

The full list is in [src/SPINGNN/KNOWN-LIMITATIONS.md](src/SPINGNN/KNOWN-LIMITATIONS.md).

## Validation

The package was compiled and checked on CPU: 185 of 185 reference dumps from
the paper's code are reproduced byte for byte, the model-call ledgers are
checked case by case, and the GPST test suites are run on three code paths
(each passes 39 of the 40 original tests; the failing test is limitation 7).
[VALIDATION.md](VALIDATION.md) lists every check, its result and what is not
covered.

## Citation

If you use this code, please cite:

```bibtex
@misc{zheng2026efficient_spin_lattice,
  author        = {Yuanmai Zheng and Liangliang Hong and Xingao Gong and Hongjun Xiang and Hongyu Yu},
  title         = {Efficient spin-lattice dynamics and spin thermometry with machine-learned magnetic potentials},
  year          = {2026},
  eprint        = {2611.32129},
  archivePrefix = {arXiv},
  primaryClass  = {physics.comp-ph},
  url           = {https://arxiv.org/abs/2611.32129}
}
```

The same reference is in [CITATION.cff](CITATION.cff), which GitHub uses for
the "Cite this repository" button. This package builds on SpinGNN; please also
cite the SpinGNN papers listed in its
[CITATION.cff](https://github.com/Hyyu-Lab/SpinGNN/blob/main/CITATION.cff).

## License

| Scope | License |
|---|---|
| `tools/gpst/` (offline GPST tools) | GNU Affero General Public License, version 3 only (`AGPL-3.0-only`); see [tools/gpst/LICENSE](tools/gpst/LICENSE) |
| Everything else: the LAMMPS package (`src/`, `patches/`, `cmake/`, `doc/`, `examples/`, `install.sh`) and the repository documentation | GNU General Public License, version 2 only (`GPL-2.0-only`); see [LICENSE](LICENSE). Code that modifies or links into LAMMPS follows the LAMMPS license. |

Original copyright and source notices are preserved.

## Questions and bug reports

Please open an issue on GitHub.
