# CVODE spray integration validation

This report preserves the initial implementation's results, including its
original order-1 default and original exchange error scales. The current default
uses BDF order 2 and exchange scales based on physical-state accuracy; see
[the tuning follow-up](TUNING.md) for the current solver's checks and timings.

The default explicit integrator preserves the baseline bit for bit in the tested
PeleLMeX spray configurations. The optional CVODE integrator completes the same
case matrix and improves the independently checked stiff Stokes problem.
An existing OpenMP grid-decomposition failure remains reproducible in both the
original and changed explicit code; it is not reported as a passing test.

## Scope and provenance

- Development checkout: `src/myBranch/PeleLMeX`; PeleC testing is deferred.
- PeleLMeX baseline: `dd28cee3e8c2c81ce9720e003655bb114b53624c`.
- PelePhysics baseline: `5e697a009eedfd9762f5dcfac60c1817d1b9971d`.
- Platform: macOS arm64, LLVM, Open MPI 5.0.11, SUNDIALS 7.2.0; 2D double
  precision, two MPI ranks. OpenMP checks additionally use two threads per rank.
- Original executables were preserved before rebuilding the modified dependency.
  Baseline builds used a separate Physics source snapshot and build caches.
  The original tracked spray sources are also archived in the evidence directory.
- GNU make builds were executed. The CMake source list was updated but a CMake
  build was not validated. GPU and 3D builds and large-particle-count performance
  are not established by these checks.

Historical validation outputs were retained locally outside the source checkout
in `src/myBranch/tests/spray_cvode`; they are not included in this repository.
The machine-readable `summary.json` contains
commands, completion status, comparison counts, grid results, and backend checks.
Every run has its own `command.json` and `run.log`. Executables and build logs
are preserved there too. `README.md` explains how to reproduce the checks.

## Default explicit regression results

The matrix covers both spray regression directories, every named SingleDropEvap
case, the upstream spray CI configurations, and the self-contained CMake-test
input. The original inputs, physical models, and case-generation utilities were
used. The comparison hashes plotfile headers and fluid FABs, ASCII particle
histories, and complete fluid/particle checkpoint state. Build metadata is
excluded because build paths and times change.

| Group | Cases/runs | Checked files | Result |
|---|---:|---:|---|
| SprayTest: two grid layouts, continuous/restart | 4 | 341 | Bitwise identical |
| Two-component GCM: WongLin, Nomura, Daif, RungeHep/Dec/Mix | 6 | 4,170 | Bitwise identical |
| 67-component GCM: RungeJP8, Burger 1/10/50 bar | 4 | 1,395 | Bitwise identical |
| One-component GCM: corresponding `-H` cases | 4 | 1,395 | Bitwise identical |
| Detailed gas mechanism: RungeJP8-D | 1 | 30 | Bitwise identical |
| Two-component MP/Antoine: six basic cases | 6 | 635 | Bitwise identical |
| Two-dimensional manifold: Daif, WongLin | 2 | 1,120 | Bitwise identical |
| Self-contained `single-drop-evap-gcm-test.inp` | 1 | 210 | Bitwise identical |
| **Total** | **28** | **9,296** | **All identical** |

Each group has a `changed/<group>/comparison.json` with empty `missing` and
`different` lists. These are comparisons against newly executed original-source
baselines, rather than experimental measurements or historical binary files.
The final rebuilt SprayTest was separately rerun and again matched all 341 files.

Run limits are explicit: SprayTest runs through step 9, with restart at step 5.
The six GCM2 cases and both default manifold cases run to their reference-data
end times: WongLin 4 s, Nomura 2.94 s, Daif 18.1 s, RungeHep 23.9 s,
RungeDec 226.6 s, RungeMix 252.1 s. Burger cases finish at 0.09 s.
JP8, JP8-H and JP8-D use the upstream JP8 CI cap of 1,000 steps (5 s).
The additional MP cases also use a 1,000-step cap; WongLin and Nomura reach
their end times first. The embedded GCM input runs to 4 s.

The original SprayTest acceptance criterion compares six variables across three
AMR levels with a maximum relative error of `1e-12`:

| MPI comparison | Explicit maximum relative error | CVODE maximum relative error |
|---|---:|---:|
| 32-box versus 64-box layout | 8.517e-15 | 4.003e-13 |
| 32-box continuous versus restart | 2.631e-15 | 2.232e-13 |
| 64-box continuous versus restart | 2.778e-16 | 1.166e-13 |

All six comparisons pass. See
`changed/final-spray/grid_checks.json`
and `cvode/final-spray/grid_checks.json`.
The CVODE grid/restart run uses default solver settings, including BDF order 1.

## Optional CVODE execution and numerical evidence

All 28 corresponding MPI runs complete with `particles.integrator=cvode_bdf`.
CVODE uses the 1,000-step cap for long GCM2, MP and manifold cases; the full
explicit runs above provide the backward-compatibility evidence. WongLin,
Nomura, Burger and the embedded input finish their requested times. GCM2 and
SprayTest use default BDF order 1; other configurations explicitly select maximum
order 2. The successful 67-component results are in `cvode/gcm67-verified`.
The final GCM2, MP, manifold and embedded-input reruns are recorded in
`cvode/final-*`. Earlier directories containing failed development attempts are retained and are
not counted as successful runs.

During validation, the Burger1Bar case exposed slightly negative component masses
in root interpolation at simultaneous depletion. The implemented event handling
now transfers those tolerance-sized residuals consistently and restarts CVODE;
all four 67-component cases subsequently complete. Accepted steps still use
CVODE's nonnegative-mass constraints.

The production RHS and integrator are also exercised against the exact solution
`u(h) = u(0) exp(-h/tau)` for a 1-micrometre droplet in the Stokes regime.
The gas is frozen, the liquid does not evaporate, and the initial speed is
`0.001 m/s`. Both methods receive the same state and requested interval.

| h/tau | Existing explicit absolute velocity error (m/s) | CVODE absolute velocity error (m/s) |
|---|---:|---:|
| 3 | 4.588e-5 | 2.572e-10 |
| 1,000 | NaN | 5.388e-23 |

These CVODE checks use `rtol=1e-9`, `atol_velocity=1e-12`, maximum order 2.
The maximum relative momentum-exchange residual is `2.726e-15`.
Saturated-gas heating satisfies `Q = m cp (Tfinal - Tinitial)` with relative
residual `1.917e-14` and reaches 350 K. Complete evaporation inside a longer
requested stage conserves mass with relative residual `6.158e-17`; this also
checks normalization by the full requested stage interval.
See the `backend-final/balances.log`.

A separate order-1 check with tighter tolerances (`rtol=1e-11`,
`atol_velocity=1e-14`, `atol_temperature=1e-8`, `max_steps=2000000`) also passes:
Stokes errors are `4.976e-10` and `5.392e-25 m/s`; the momentum and heating
residuals are below `9e-14` and `3e-13`, respectively. Its command is in
`README.md`, and the log is `relaxation-final-order1.log` in the evidence directory.
The larger step budget is for this tight reference check, not the production
default.

Unknown integration methods, negative tolerances, and an exhausted one-step
budget all produce the expected clear failures. The
`backend-final/backend_checks.json`
contains commands and return codes. The improvement established here is stiff
ODE stability and accuracy with consistent exchange accounting. Improved
agreement with every experimental evaporation dataset is not claimed.

## OpenMP limitation reproduced in the baseline

The original `SprayParticles.H` OpenMP pragma uses unqualified
`Gpu::notInLaunchRegion()` outside the AMReX namespace and did not compile with
OpenMP. Qualifying it as `amrex::Gpu::notInLaunchRegion()` was necessary for both
the baseline and changed OpenMP builds. This is the only source adjustment made
to the original OpenMP baseline.

With two MPI ranks and two OpenMP threads, the original and changed explicit
executables produce **341 bitwise-identical files**, including checkpoints.
Both fail the original 32-box/64-box grid comparison by exactly the same amount:
maximum relative error `0.09930719373`. Both restart comparisons pass, with maxima
`2.844e-16` and `2.778e-16`.

The CVODE OpenMP execution also completes and its restart checks pass
(`3.067e-13` and `3.108e-13`); its layout comparison fails (`0.005420331103`).
This is not evidence that the inherited grid dependence is resolved. Investigation
of the existing OpenMP particle/deposition behavior is separate work.

See the `openmp/explicit-spray/comparison.json`,
`openmp/baseline-explicit-spray/grid_checks.json`,
`openmp/explicit-spray/grid_checks.json`,
and `openmp/cvode-spray/grid_checks.json`.
