# Exchange error scales and BDF-order follow-up

This follow-up changes only the optional CVODE particle integrator. The production
default remains `particles.integrator=explicit`. Development and validation use
`src/myBranch/PeleLMeX`; PeleC testing remains deferred.

## Changes

The default maximum BDF order is now 2, allowing CVODE to increase its order after
startup. `particles.cvode.max_order=1` still selects backward Euler; values 3–5
remain available. Particle-state tolerances retain their original defaults.

The exchange accumulators still start at zero and participate in CVODE's error
test, but their absolute tolerances now include physical reference scales:

```text
dm       = Nfuel * atol_mass + rtol * initial_mass
u_ref[j] = max(abs(initial_particle_velocity[j]), abs(gas_velocity[j]))
du[j]    = atol_velocity + rtol * u_ref[j]
atol_J[j] = initial_mass * du[j] + dm * u_ref[j]
```

The analogous energy scale uses liquid heat capacity and
`atol_temperature + rtol * max(initial_particle_temperature, gas_temperature)`.
Its mass-error contribution uses the larger sensible-energy/vapor-enthalpy
reference. The PeleC branch additionally propagates kinetic-energy uncertainty;
that branch has not been tested here.

The source equations, constraints, depletion handling, particle-state updates,
and deposition are unchanged. Each particle/stage still reinitializes the reused
workspace: numerical history cannot be shared between different particles.

## Correctness checks

Validation outputs were retained locally under `src/myBranch/tests/spray_cvode/tuning`
and are not included in this repository.

- All **28 revised-default CVODE MPI runs complete** across SprayTest, GCM2,
  GCM67, GCM1, detailed gas mechanism, MP2, manifold2, and the embedded GCM input.
  Named SingleDropEvap cases use a 1,000-step cap; short cases reach their end time
  first. The embedded input runs to its original end time.
  Commands and completion counts are in `regression-summary.json`.
- The original SprayTest grid/restart criterion remains **1e-12**. Maximum
  relative errors are **9.697e-13** for the layout comparison, **3.055e-13** for
  the 32-box restart, and **2.993e-13** for the 64-box restart. See
  `tuning/cvode/spray/grid_checks.json`.
- The updated executable's default explicit SprayTest rerun produces **341
  bitwise-identical files**, including checkpoints. See
  `tuning/explicit/spray/comparison.json`.
  The complete original explicit matrix is preserved in `VALIDATION.md`.
- Exact-solution, heating, evaporation, and expected-failure checks all pass in
  `tuning/backend-verified/backend_checks.json`.

The exact-solution driver now includes a non-evaporating **1 mm droplet at
15 m/s**, matching SprayTest's initial droplet size/speed. In the intermediate
Reynolds regime, the drag equation has the form
`du/dt = -a u (1 + b u^(2/3))`, with a closed-form solution. The measured absolute
velocity error is **9.097e-8 m/s**, and the relative momentum-balance residual is
**5.107e-16**.

These reference checks use `rtol=1e-10`, `atol_velocity=1e-12`,
`atol_temperature=1e-8`, and maximum order 2. The initial nonlinear check at the
previous reference tolerance, `rtol=1e-9`, produced `1.052e-6 m/s`, slightly above
its fixed `1e-6 m/s` bound. The reference tolerance was tightened; the accuracy
bound was not weakened. The initial failure is retained in `tuning/backend`.

Other final reference results:

| Check | Result |
|---|---:|
| Stokes velocity error, h/tau=3 | 1.789e-10 m/s |
| Stokes velocity error, h/tau=1,000 | 1.927e-23 m/s |
| Maximum Stokes momentum-balance residual | 6.902e-15 |
| Saturated-gas heating energy residual | 4.219e-14 |
| Heating final temperature | 350.0000004 K |
| Complete-evaporation mass residual | 6.158e-17 |

## Serial repeated timing comparison

`benchmark.py` uses two MPI ranks and the same first flow step, with physical
interval **0.005208333333 s**, for every variant. Three repetitions run serially,
rotating variant order. Other agent-run builds and simulations finish before
these measurements. Executable hashes and every command are retained. Times are
wall-clock medians, and advance time excludes initialization and final output.

For the 32-box layout:

| Exchange scales | Maximum BDF order | Advance time | Whole-program time |
|---|---:|---:|---:|
| Original | 1 | 17.485447 s | 17.837494 s |
| Original | 2 | 1.522143 s | 1.885660 s |
| Physical state | 1 | 2.756024 s | 3.116526 s |
| Physical state | 2, new default | 0.871270 s | 1.234458 s |
| Native explicit | — | 0.065941 s | 0.430777 s |

The combined change reduces advance time by **20.1 times**, and whole-program
time by **14.4 times**, relative to the previous implicit default. At order 1,
scaling alone reduces advance time by **6.34 times**; at order 2 its improvement
is **1.75 times**. Higher order with the original scales improves advance time
by **11.5 times**. The effects overlap, so these factors should not be multiplied.
The revised implicit advance is still **13.2 times slower than explicit** for
this mildly stiff initial step. Implicit integration has not been made a
universal speed improvement.

See `tuning/benchmark32/timings.json`
and the accompanying per-run logs. These results apply to this CPU, mesh layout,
state, and tolerances; they are not cross-platform performance guarantees.

The independently repeated 64-box layout gives:

| Exchange scales | Maximum BDF order | Advance time | Whole-program time |
|---|---:|---:|---:|
| Original | 1 | 32.218632 s | 32.593438 s |
| Original | 2 | 2.736467 s | 3.108619 s |
| Physical state | 1 | 5.129650 s | 5.497904 s |
| Physical state | 2, new default | 1.512192 s | 1.881069 s |
| Native explicit | — | 0.095345 s | 0.468249 s |

That is **21.3 times faster** in advance time and **17.3 times faster** in
whole-program time than the original implicit default. The revised advance is
still **15.9 times slower than explicit**. See
`tuning/benchmark64/timings.json`.
Every run in both layouts finishes at the same physical time.

## Accuracy against tighter coupled-flow references

The 32-box first step was also run with maximum BDF order 5 and three successively
tighter tolerance sets. The final reference uses `rtol=1e-11`,
`atol_temperature=1e-10`, `atol_velocity=1e-12`, and `atol_mass=1e-26`.
Reference commands are retained in `reference9`, `reference10`, and `reference11`.
All 65,536 particles were matched by particle ID and originating MPI rank.

Maximum absolute particle errors relative to `reference11`:

| Variant | Velocity error (m/s) | Temperature error (K) | Diameter error (m) |
|---|---:|---:|---:|
| Explicit | 2.261e-2 | 2.054e-2 | 1.477e-9 |
| Original scales, order 1 | 5.045e-5 | 4.689e-5 | 3.515e-12 |
| Original scales, order 2 | 2.906e-6 | 1.008e-6 | 1.546e-13 |
| Physical scales, order 1 | 3.393e-4 | 3.145e-4 | 2.346e-11 |
| Physical scales, order 2, new default | 4.317e-5 | 1.429e-5 | 3.097e-12 |
| Reference10 versus reference11 | 8.776e-8 | 2.855e-8 | 4.872e-15 |

Scaling alone at order 1 trades accuracy for speed. The combined default restores
and slightly improves accuracy relative to the previous order-1 default for this
step, while reducing its advance time by 20.1 times. The original order-2 setting
remains more accurate than the revised default in this comparison; it also
requests stricter exchange accuracy and costs more. The tighter-reference
uncertainty is small compared with the reported default errors.

Gas fields include the effect of deposited particle sources. Maximum relative
errors over all three AMR levels are:

| Field | Explicit | Previous implicit default | Revised implicit default | Reference10 versus reference11 |
|---|---:|---:|---:|---:|
| Fuel partial density | 3.645e-2 | 1.054e-4 | 9.414e-5 | 1.842e-7 |
| Gas x velocity | 1.580e-2 | 4.481e-5 | 2.933e-5 | 9.449e-8 |
| Gas temperature | 1.598e-7 | 9.366e-10 | 4.283e-10 | 1.110e-12 |
| Enthalpy density | 3.441e-7 | 2.008e-9 | 8.114e-10 | 3.547e-12 |

See `tuning/particle-accuracy.json`
and `tuning/gas-accuracy/comparisons.json`.
These are solver-tolerance convergence checks at fixed spatial resolution and
flow timestep, rather than a spatial/flow-timestep convergence study or a claim
of improved agreement with every experimental dataset.

The combined machine-readable record, including source/executable hashes, is
`tuning/summary.json`.

The inherited OpenMP layout-dependence described in `VALIDATION.md` is outside
this change. GPU, 3D, CMake, and PeleC validation remain outside this follow-up.

## Explicit-accuracy comparison

For the same one-step, two-rank, 65,536-particle case at `t=0.005208333333 s`, the table shows maximum particle errors against `reference11` and median timings over three runs per variant. Both BDF trials use `rtol=0.1` and the compiled absolute tolerances; production defaults remain unchanged. Relative mass error uses diameter cubed.

| Integrator | Max velocity error (m/s) | Max temperature error (K) | Max diameter error (m) | Max relative mass error | Advance (s) | Total (s) |
|---|---:|---:|---:|---:|---:|---:|
| Explicit | 0.0226067 | 0.0205416 | 1.47688e-9 | 4.43075e-6 | 0.064565 | 0.430036 |
| BDF, max order 1 | 0.0114100 | 0.0107490 | 8.21424e-10 | 2.46434e-6 | 0.431121 | 0.801517 |
| BDF, max order 2 | 0.00401016 | 0.00377224 | 2.87214e-10 | 8.61663e-7 | 0.446386 | 0.825326 |

The order-2 velocity error was 0.00523243 m/s at `rtol=1e-3`, and 0.00401016 m/s at both `4e-3` and `0.1`; loosening tolerance alone did not reach the explicit error in these trials. Order 1 was tested to obtain an approximate match: its velocity error is within a factor of two of explicit, though still smaller. Its median advance time is 6.68 times the explicit time for this case. These measurements do not establish a general accuracy floor or performance ratio.

Maximum `fcompare` relative errors across the three levels, using each candidate as the denominator, are below. These are field norm ratios, not cellwise relative maxima; first runs were compared against `reference11`.

| Integrator | Fuel partial density | Gas x velocity | Temperature | rhoh | Density | Spray volume |
|---|---:|---:|---:|---:|---:|---:|
| Explicit | 3.64548e-2 | 1.58021e-2 | 1.59751e-7 | 3.44111e-7 | 8.87973e-9 | 3.96267e-6 |
| BDF, max order 1 | 1.92895e-2 | 8.21103e-3 | 8.52276e-8 | 1.82530e-7 | 4.91866e-9 | 2.20303e-6 |
| BDF, max order 2 | 6.81869e-3 | 2.87420e-3 | 2.98216e-8 | 6.38803e-8 | 1.71970e-9 | 7.70280e-7 |
