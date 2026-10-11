# Running the Three Cases

English | [简体中文](README.zh-CN.md)

## 1. Cases and Inputs

| Directory | Mesh cells | Configurations | Time step / s | Full steps | Write interval / steps |
| --- | ---: | --- | ---: | ---: | ---: |
| `couette/` | 500 | Fixed weight, MCOR exact collisions, MCOR stochastic updates | `1e-5` | 40000 | 1000 |
| `plate/` | 6000 | Fixed weight, MCOR exact collisions | `4e-6` | 40000 | 100 |
| `orifice/` | 135168 | Fixed weight, MCOR exact collisions | `1e-7` | 4000 | 100 |

There are 7 input configurations. Each directory contains `constant`, `system`, and `case.foam`; meshes are generated from dictionaries. The locally refined orifice configuration has a 2 mm aperture diameter, a 12 mm reservoir diameter, an axial length of 3 mm on each side, and a pressure ratio of 100:1.

Each model configuration has one standard input directory. Directory names use `fixed` for the fixed-weight control, `moment_exact` for MCOR exact collisions, and `moment_stochastic` for MCOR stochastic updates. Use `--ranks` to set the process count at run time.

The inputs use the solver's default random-number initialization, based on the current time and a process-dependent offset. Assess statistical variability through independent runs in separate output directories.

## 2. Environment and Run Command

Complete the build using the [installation guide](../README.md), then load the environment in Bash and create the parent output directory:

```bash
source /opt/OpenFOAM/OpenFOAM-v1706/etc/bashrc
cd "$HOME/OpenFOAM/MCOR1_0"
source etc/bashrc
mkdir -p "$HOME/MCOR_runs"
```

The run command is:

```text
./cases/Allrun FAMILY/CASE NEW_OUTPUT [--ranks N] [--steps N] [--prepare-only]
```

- The script creates the output directory. Choose a new name for each run.
- The default is a serial run using the full duration and write interval from the input files.
- `--ranks N` sets the MPI process count and updates the decomposition dictionary and parallel-regulation option in the output copy. MPI inherits process-binding settings from the current environment.
- `--steps N` selects a short test with the original time step, N total steps, and one statistical-field write at the final step.
- `--prepare-only` performs mesh generation, initialization, and any required decomposition, then finishes with the initialized case in the specified output directory.

Required tools are OpenFOAM's `blockMesh`, `checkMesh`, and `foamDictionary`, plus this package's `dsmcInitialise+` and `dsmcFoam+`. The orifice case also uses `topoSet`, `refineMesh`, and `createBaffles`. Parallel runs require `mpirun`, `decomposePar`, `reconstructPar`, and Scotch decomposition support.

## 3. Short Tests

The following commands run the MCOR branch of each case. Adjust the OpenFOAM and package paths to your installation, then execute from the package root:

```bash
./cases/Allrun couette/couette_moment_exact "$HOME/MCOR_runs/couette_check" --steps 10
./cases/Allrun plate/plate_moment_exact "$HOME/MCOR_runs/plate_check" --ranks 2 --steps 5
./cases/Allrun orifice/orifice_moment_exact "$HOME/MCOR_runs/orifice_check" --ranks 2 --steps 1
```

For the fixed-weight branch, replace `moment_exact` with `fixed` and choose a new output directory. The short orifice test uses the full 135168-cell mesh and initializes approximately several million simulation parcels. Check available memory and disk space before running it.

Successful execution ends with `CASE_RUN_COMPLETE`. Inspect `log.checkMesh`, `log.dsmcInitialise`, and `log.dsmcFoam` in the output directory. MPI runs also produce `log.decomposePar` and `log.reconstructPar`. See the [test guide](../tests/README.md) for coverage.

## 4. Full Runs

Allocate CPU, memory, and disk resources according to your computing platform, then set `--ranks` to the allocated process count. These commands retain the full input durations and original process counts: 1 for Couette, 60 for plate, and 32 for orifice.

```bash
./cases/Allrun couette/couette_fixed "$HOME/MCOR_runs/couette_fixed" --ranks 1
./cases/Allrun couette/couette_moment_exact "$HOME/MCOR_runs/couette_mcor" --ranks 1

./cases/Allrun plate/plate_fixed "$HOME/MCOR_runs/plate_fixed" --ranks 60
./cases/Allrun plate/plate_moment_exact "$HOME/MCOR_runs/plate_mcor" --ranks 60

./cases/Allrun orifice/orifice_fixed "$HOME/MCOR_runs/orifice_fixed" --ranks 32
./cases/Allrun orifice/orifice_moment_exact "$HOME/MCOR_runs/orifice_mcor" --ranks 32
```

For the Couette stochastic-update branch, select `couette/couette_moment_stochastic` and choose a new output directory.

## 5. Execution and Outputs

`Allrun` copies the selected inputs into a new output directory, then performs these steps:

1. Checks this package's applications and dynamic libraries, sets the process count and optional short-test duration, and writes `run-settings.txt`.
2. Generates the mesh. The orifice workflow applies outer-region refinement, near-aperture refinement, region selection, and baffle creation in that order.
3. Runs `checkMesh` and checks for `Mesh OK.`.
4. Initializes parcels serially and writes the collision-kernel upper bound from `constant/cpcSigmaTcRMax` into `0/dsmcSigmaTcRMax`.
5. Runs `decomposePar` for parallel cases, then starts the solver with the requested process count.
6. Checks for `End main` and reconstructs the final time after MPI runs.

Parallel sampling times are stored in `processor*` directories. To reconstruct all sampled fields, run:

```bash
reconstructPar -case "$HOME/MCOR_runs/plate_mcor" -noLagrangian -noZero
```

Open `case.foam` from the output directory in ParaView to inspect the mesh and statistical fields. For quantitative comparisons, use the original sampling windows: steps 30001--40000 for Couette, 6001--40000 for plate, and 2001--4000 for orifice.

## 6. Calculation Settings

All 7 input configurations use MCOR1_0. Keep the mesh, time step, gas conditions, MPI decomposition, and sampling procedure consistent when comparing methods.

Research calculations require sampling over the stated durations and windows, together with the necessary statistical-accuracy, discretization, and physical-model checks.
