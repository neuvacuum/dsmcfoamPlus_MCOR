# MCOR1_0

English | [简体中文](README.zh-CN.md)

MCOR1_0 builds and runs the dsmcFoam+ moment-conservative occupancy regulation (MCOR) extension on **OpenCFD OpenFOAM v1706**.

The package includes moment-conservative occupancy regulation, weighted collisions and statistics, reference-weight free-stream injection, and Couette, plate, and orifice cases. See the [implementation guide](docs/implementation.md) for configuration and [NOTICE](NOTICE.md) for licensing and third-party notices.

## 1. Requirements

- Linux and Bash.
- An installed OpenCFD OpenFOAM v1706 environment with `wmake`, source headers, compiled base libraries, and `blockMesh`.
- GNU make, a compatible C++ compiler, MPI, and third-party libraries with matching ABI settings.
- Standard Linux utilities including `ldd`, `awk`, `grep`, `tee`, and `mktemp`.
- A user-writable installation directory, preferably using letters, digits, and underscores, such as `$HOME/OpenFOAM/MCOR1_0`.

The tested platform is `linux64GccDPInt32Opt`: double precision, 32-bit labels, GCC, and Open MPI. See the [test guide](tests/README.md) for test coverage and commands.

## 2. Installation

Run the following commands as a regular user in Bash. Place the complete package folder at `$HOME/OpenFOAM/MCOR1_0`, or adjust the paths to your installation location.

### 2.1 Load OpenFOAM

Open a new terminal and replace the path in the first command with your OpenFOAM installation path:

```bash
source /opt/OpenFOAM/OpenFOAM-v1706/etc/bashrc
echo "$WM_PROJECT_VERSION"
echo "$WM_OPTIONS"
command -v wmake
command -v blockMesh
```

Confirm that the reported OpenFOAM release is `v1706` and that both tools belong to this installation.

### 2.2 Load and Build the Extension

```bash
cd "$HOME/OpenFOAM/MCOR1_0"
source etc/bashrc
./Allwmake 2
```

The argument `2` sets the number of compilation jobs. Use `./Allwmake 1` when memory is limited. The script builds three libraries and two applications, saving logs to `logs/build.*.log`. A successful build ends with `MCOR_INSTALLATION_PASS` and `MCOR_BUILD_PASS`.

Build products are stored inside the package:

```text
platforms/linux64GccDPInt32Opt/bin/dsmcFoam+
platforms/linux64GccDPInt32Opt/bin/dsmcInitialise+
platforms/linux64GccDPInt32Opt/lib/liblagrangian+.so
platforms/linux64GccDPInt32Opt/lib/libgeneralMolecule.so
platforms/linux64GccDPInt32Opt/lib/libdsmcFoam+.so
```

The actual platform directory is determined by `$WM_OPTIONS`. The package's `etc/bashrc` points `WM_PROJECT_USER_DIR`, `FOAM_USER_APPBIN`, and `FOAM_USER_LIBBIN` to this package and prepends its paths to `PATH` and `LD_LIBRARY_PATH`. Use a separate terminal and the corresponding environment for other solver installations or user extensions.

### 2.3 Check the Installation and Run a Smoke Test

```bash
./tests/check-installation.sh
./tests/run-smoke.sh
```

The first command checks that both applications start and that all three extension libraries load from this package. The second generates a mesh, initializes the gas, and runs a short case in a new `runs/smoke.*` directory. Success is reported as `MCOR_SMOKE_PASS`, followed by the output directory.

The smoke test checks 10 time steps, 20 regulation-stage records, nonzero candidate and accepted collisions, and pruning operations. The maximum scaled moment residual in the log must be at most `1e-9`. A summary is written to `smoke-summary.txt` in the output directory.

To specify a new output directory for the script to create:

```bash
./tests/run-smoke.sh "$HOME/mcor-smoke-01"
```

These checks verify installation, dynamic linking, and basic execution of the supplied serial case. Research calculations require case-specific mesh and time-step convergence checks, statistical accuracy assessment, and physical-model validation. Parallel calculations also require MPI checks.

## 3. Usage

For each new terminal, load the base OpenFOAM environment followed by this package:

```bash
source /opt/OpenFOAM/OpenFOAM-v1706/etc/bashrc
source "$HOME/OpenFOAM/MCOR1_0/etc/bashrc"
command -v dsmcFoam+
command -v dsmcInitialise+
```

Back up your case inputs and results before working on an existing case. Use `examples/closedBox` as a reference for the required dictionaries, and set gas density, time step, and mesh for your physical problem.

For a new case, run:

```bash
blockMesh -case /absolute/path/to/new-case
dsmcInitialise+ -case /absolute/path/to/new-case
dsmcFoam+ -case /absolute/path/to/new-case
```

Check boundary models, species, reference parcel weights, and `adaptiveProperties` against the [implementation guide](docs/implementation.md). For parallel execution, prepare a compatible MPI environment and decomposition dictionary, set `allowRuntimeMomentParallel true`, and perform the corresponding parallel tests.

The `cases/` directory provides 7 input configurations across the Couette, plate, and orifice cases. Follow the [case guide](cases/README.md) to use `cases/Allrun` for mesh generation, initialization, serial or MPI solving, and final-time reconstruction.

## 4. Troubleshooting

| Symptom | Action |
| --- | --- |
| A message requests the v1706 environment | In a new terminal, source the correct OpenFOAM `etc/bashrc`, then this package's environment |
| Missing headers or `wmake` | Complete the OpenFOAM v1706 development installation |
| `libdsmcFoam+.so` cannot be found | Source this package's `etc/bashrc` and run `tests/check-installation.sh` to check the resolved paths |
| `undefined symbol` or `GLIBCXX` errors | Check compiler and OpenFOAM/MPI ABI compatibility and library paths; rebuild in a compatible environment |
| A newer compiler treats legacy-interface warnings as errors | Use the tested build environment and inspect the compiler diagnostics |
| Compilation is killed by the system | Use `./Allwmake 1` and check available memory |
| The dynamic collision branch exits because acceptance exceeds 1 | Correct the initial collision bound and case settings, then rerun |
| Particle counts differ from their targets | Check the occupancy band and assess pruning feasibility and internal-state constraints |

After changing the compiler, precision, label width, or OpenFOAM build, compile a fresh source copy. Close the terminal when finished using the environment. Back up the required results in `runs` before moving or removing the package.

## 5. Publishing to GitHub

Use the contents of `MCOR1_0/` as the GitHub repository root, including the hidden `.gitignore` file. The repository homepage will display the English `README.md`, with a link to the Chinese translation at the top.

Commit the source code, all three cases, installation scripts, documentation, `VERSION`, `COPYING`, and `NOTICE.md` together. Attach `MCOR1_0.tar.gz` to a GitHub Release as the downloadable source archive.
