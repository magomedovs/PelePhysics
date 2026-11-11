# OpenMP layout discrepancy: confirmed cause

## Finding

The observed layout-dependent spray solution is caused by premature finalization
of the shared particle-timestep reduction in
`PelePhysics/Source/Spray/SprayParticles.cpp`, in `estTimestep()`.

Particle contributions are evaluated in an OpenMP region. The original routine
also calls `reduce_data.value()` and updates shared `dt` inside that region,
after each thread finishes its own iterator. There is no guarantee that all
other threads have finished contributing at that point.

The CPU AMReX reduction has per-thread partial-result slots. Its `value()` merges
these slots into the first slot and sets `m_result_is_ready`. Subsequent `value()`
calls return the cached result. A worker can therefore finalize a partial minimum
before another worker contributes a smaller particle timestep. The later
contribution is not incorporated into the cached global result. Concurrent
finalization and writes to `dt` also constitute data races.

Missing the true minimum permits a timestep that is too large for the fastest
parcel. It can change the motion subcycles, particle trajectories, and gas
exchange. Changing the box layout or regridding changes the division and timing
of work between threads, exposing the defect. This is separate from CVODE's
internal adaptive stepping.

Relevant source:

- [PelePhysics particle timestep routine](../../Source/Spray/SprayParticles.cpp),
  especially the parallel region and its final reduction near lines 104–132.
- [AMReX CPU reduction implementation](https://github.com/AMReX-Codes/amrex/blob/5a9a0cc7f8d8/Src/Base/AMReX_Reduce.H),
  `ReduceData`, thread-local `reference(tid)`, and cached `ReduceOps::value()`.

## Isolated correction

Keep particle-rate evaluation parallel. Move reduction finalization and the
single update to `dt` after the end of the parallel region:

```cpp
#pragma omp parallel if (Gpu::notInLaunchRegion())
{
  for (MyParConstIter pti(*this, level); pti.isValid(); ++pti) {
    // Existing reduce_op.eval(...) calls, unchanged.
  }
}
ReduceTuple hv = reduce_data.value();
Real ldt_cpu = amrex::get<0>(hv);
dt = amrex::min<Real>(dt, ldt_cpu);
```

The end of the parallel region synchronizes the workers before the reduction is
read. Jet-velocity checks and the MPI minimum reduction remain unchanged.
Protecting only the assignment to `dt` would not address premature cached
finalization.

The correction was tested in an isolated Physics source overlay. The production
working-tree routine has **not** been patched by this investigation. The exact
`omp_investigation/reduction-fix.patch`
is retained for review/application. The pre-investigation executable was restored
after testing, and the diagnostic cached objects were invalidated.

## Controlled evidence

All runs use the existing 2D SprayTest case, two MPI ranks, and the original
32-box/64-box layouts and checkpoint/restart checks. No tolerances, physics, or
spatial resolutions were changed. The isolated correction moves only the three
final-reduction statements.

Before the correction, the retained original and changed explicit OpenMP runs
have the same failure. At two threads per rank, the largest AMR-level relative
field-norm discrepancy is **0.09930719373** in `spray_vol`. This is a comparison
of spatial fields, not a measurement of total mass loss.

The first logged timestep discrepancy occurs at step 4, around regridding:

| Configuration | 32-box dt | 64-box dt | Final time after nine steps |
|---|---:|---:|---|
| Original routine, two threads | 0.006224799686 | 0.006222598231 | 0.05657581533 versus 0.05654482644 |
| Corrected routine, two threads | 0.006222598231 | 0.006222598231 | Both 0.05654482644 |

The initial state and logged timesteps through the first four advances agree;
the correction restores the same step-4 timestep and later final time as MPI-only.

Each test group below contains four completed runs: continuous/restart for both
layouts. All comparison errors are maximum relative field-norm errors over the
six original test variables and three AMR levels.

| Routine/integrator | Threads per rank | Layout error | 32-box restart | 64-box restart |
|---|---:|---:|---:|---:|
| Original explicit, control | 1 | 8.517e-15 | 2.631e-15 | 2.778e-16 |
| Corrected explicit | 2 | 8.517e-15 | 2.631e-15 | 2.778e-16 |
| Corrected explicit | 4 | 8.517e-15 | 2.631e-15 | 2.778e-16 |
| Corrected current CVODE default | 2 | 9.697e-13 | 3.055e-13 | 2.993e-13 |

Every corrected check passes the original **1e-12** criterion. Additionally,
**all 341 checked files from the corrected two-thread explicit run match the
MPI-only original-source baseline bit for bit**, including fluid and particle
checkpoint state. Build metadata is excluded. The CVODE comparison results equal
those of the current MPI-only CVODE validation.

This provides a controlled explanation and correction for the observed failure.
It is not a claim that every possible OpenMP configuration is free of other bugs.
GPU builds and PeleC were not tested in this investigation.

## Reproduction artifacts

`omp_investigation/summary.json`
contains the complete check results. The same directory retains:

- `reduction-fix.patch` and the diagnostic `SprayParticles` source;
- the pre-correction and corrected executables and the build log;
- per-run commands, logs, checkpoints, and plotfiles;
- `mpi_comparison.json`, listing the 341-file bitwise comparison;
- the original one-thread control and corrected 2-/4-thread and CVODE suites.

The original failing results remain in
`openmp/explicit-spray/grid_checks.json`.
The overlay build uses GNU make with `COMP=llvm USE_MPI=TRUE USE_OMP=TRUE` and
`OMPI_CXX=/opt/homebrew/opt/llvm/bin/clang++`; the build log records its exact
compilation and link commands. Each regression run records `OMP` thread count in
its log and the MPI launch/input overrides in `command.json`.

The correction remains separate from the CVODE implementation packaged here;
this note records the isolated downstream investigation and its local artifacts.
