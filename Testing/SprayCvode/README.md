# Spray CVODE downstream validation

The [implementation report](../../Docs/SprayCvodeImplementation.md) describes the
architecture, limitations, and measured speed. The default particle integrator
remains `explicit`; the optional CPU integrator is `cvode_bdf`, with maximum BDF
order 2. Runtime controls are documented in [Spray.rst](../../Docs/sphinx/Spray.rst).

These tools exercise PelePhysics through **PeleLMeX's spray regression cases**.
They are not standalone PelePhysics executables. They use the downstream case
inputs and `CaseInfo.py`, without changing the downstream source checkout.

## Prerequisites and build

Set `PELELMEX_HOME` to a compatible PeleLMeX checkout. If it is unset, the tools
look for a sibling checkout named `PeleLMeX`. The recorded tests used PeleLMeX
`dd28cee3e8c2c81ce9720e003655bb114b53624c`, 2D double precision, GNU make,
LLVM, Open MPI 5.0.11, and SUNDIALS 7.2.0 on macOS arm64.

Run from the PelePhysics root. Build the downstream executable against this
dependency, using the compiler and installed third-party libraries appropriate
to your machine:

```sh
export PELELMEX_HOME=/path/to/PeleLMeX
export PELE_PHYSICS_HOME="$PWD"
make -C "$PELELMEX_HOME/Exec/RegTests/SprayTest" -j8 \
  COMP=llvm USE_MPI=TRUE PELE_PHYSICS_HOME="$PELE_PHYSICS_HOME" \
  TMP_BUILD_DIR=tmp_build_spray_cvode KEEP_BUILDINFO_CPP=TRUE
```

Initialize/build the dependency's required AMReX/SUNDIALS components beforehand
as described by the downstream build. Use separate caches and preserve original
executables for before/after comparisons. Serialize builds sharing a case
directory because they generate `AMReX_buildInfo.cpp` there. CMake consumers with
explicit spray source lists must add `SprayCvode.cpp` (and the two new headers).

## Regression and numerical checks

```sh
python3 Testing/SprayCvode/regression.py --group spray --np 2 \
  --exe /path/to/spray.ex --output /path/to/explicit
python3 Testing/SprayCvode/regression.py --group spray --np 2 \
  --exe /path/to/spray.ex --output /path/to/implicit --integrator cvode_bdf
python3 Testing/SprayCvode/regression.py --output /path/to/explicit \
  --compare /path/to/original-explicit
python3 Testing/SprayCvode/regression.py --output /path/to/implicit \
  --check-grid /path/to/fcompare.ex
```

Use fresh output directories. Comparisons hash plot/particle data and complete
fluid/particle checkpoint state; build metadata is excluded. Grid/restart checks
apply the upstream six-variable relative-error limit of `1e-12`.

SingleDropEvap groups require separately compiled executables:

| Group | Compile-time configuration |
|---|---|
| `gcm2` / `input-gcm` | `SPRAY_GCM=TRUE SPRAY_FUEL_NUM=2` |
| `gcm67` | `SPRAY_GCM=TRUE SPRAY_FUEL_NUM=67` |
| `gcm1` | `SPRAY_GCM=TRUE SPRAY_FUEL_NUM=1` |
| `detailed67` | GCM67 plus `Chemistry_Model=fuellib_posf_nonreacting` |
| `mp2` | `SPRAY_GCM=FALSE SPRAY_FUEL_NUM=2` |
| `manifold2` | GCM2 plus `USE_MANIFOLD=TRUE Manifold_Dim=2` |

For named long cases, `--max-steps 1000` matches the recorded CI-style limit.
Manifold generation additionally requires `--cmlm /path/to/cmlm` and CMLM's Python
dependencies; recorded tests used CMLM 0.0.1 and deterministic component ordering.
The other harness operations require only the Python standard library.

The exact-solution driver uses the production RHS and integrator for Stokes and
nonlinear drag, momentum/energy balances, and complete evaporation. Link it using
a completed **SprayTest** GNU make log containing compilation and link commands:

```sh
python3 Testing/SprayCvode/build_relaxation.py /path/to/build.log \
  /path/to/relaxation.ex
python3 Testing/SprayCvode/regression.py --output /path/to/backend-checks \
  --check-backend /path/to/relaxation.ex
```

The helper currently expects `mpicxx` commands in the log. The reference checks
use `rtol=1e-10`, `atol_velocity=1e-12`, `atol_temperature=1e-8`, and maximum
order 2. Unknown-method, invalid-tolerance, and exhausted-step-budget checks are
also included.

## Timing comparison

```sh
python3 Testing/SprayCvode/benchmark.py --before /path/to/old-spray.ex \
  --after /path/to/new-spray.ex --output /path/to/benchmark \
  --repeats 3 --steps 1 --box-size 32 --np 2
```

The before executable must use the original exchange scales. The benchmark runs
variants sequentially, rotating their order, and records commands, hashes, median
advance/total times, and physical end times. One step gives an identical physical
interval; longer coupled runs can select different flow timesteps. Finish other
builds/simulations before measuring timings.

[VALIDATION.md](VALIDATION.md) preserves the initial compatibility results and
known OpenMP issue. [TUNING.md](TUNING.md) records the revised default and the
approximate explicit-accuracy comparison. Their referenced measurement artifacts
were retained locally outside this repository and are not distributed here.
