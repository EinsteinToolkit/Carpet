# Ticket 2920

https://github.com/EinsteinToolkit/tickets/issues/2920

This branch is `brandt/ticket-2920`. It contains `1bca6d19` ("Carpet: size subcycling buffers for the time refinement ratio") and the changes below. `CarpetLib/src/ggf.cc` still has a local time-interpolation print. That print is not on this branch. The compiled `buffer_zones_test` binary is not on this branch either. `bugfix/ticket-2920` stays at `1bca6d19`.

## What changed here

`1bca6d19` sizes a subcycling buffer from the time-refinement ratio. The width is `time_substeps * (integrator_substeps * ghosts + additional_buffer_zones) - ghosts`. With `use_tapered_grids` left off, `time_substeps` is 1, so one RK4 step with 3 ghosts is 9 zones. With it on, `time_substeps` is this level's time factor divided by its parent's, so a 2:1 RK4 level is 21 zones. A level that shares its parent's timestep stays one step wide. Adams-Bashforth 3 with one integrator substep and tapered grids off has width 0. Ghost size stays 3.

`use_tapered_grids` only chooses that width. The default is `no`. It does not choose whether a fine step refills its ghost and buffer zones.

Two parameters choose the refill:

- `use_time_interpolation`, default `yes`. On a fine step that falls between parent times, refill the ghost and buffer zones by time interpolation.
- `use_final_copy`, default `yes`. On the fine step that lands on the parent time, refill those zones from the parent at that time.

Both defaults on is the proposed schedule: every fine step prolongates, the halfway step interpolates, and the parent-time step copies. Turning either refill off makes `SetTaperForSubcycling` skip that prolongation. If the parent stride is not a multiple of this level's stride, that path aborts. A level with the same timestep as its parent has only the copy step.

`poison_taper`, default `no`, writes NaN into the halo shell an integrator substep consumed when that step did not prolongate. A later boundary prolongation clears the shell counter. The schedule entry is inside `if (poison_taper)`.

The test `buffer_subcycling_legacy` covers 2:1 RK4 with `use_tapered_grids = no`, where the width stays 9. `buffer_subcycling_evolve` keeps the 21-zone grid, skips the halfway interpolation, and copies on the parent-time step. Its `phi.x.asc` was regenerated for that copy.

## The other branches

Check out `brandt/ticket-2920` in each of these. The thorn list `einsteintoolkit.th` in the ticket 2920 directory does that for these five and leaves the rest of the Einstein Toolkit on `ET_2026_05`.

- Carpet: https://github.com/einsteintoolkit/Carpet.git
- cactusnumerical: https://bitbucket.org/cactuscode/cactusnumerical.git
- einsteinbase: https://bitbucket.org/einsteintoolkit/einsteinbase.git
- lean_public: https://bitbucket.org/canuda/lean_public.git
- Proca: https://bitbucket.org/canuda/Proca.git

Charged black hole parameter files in that same ticket directory: `charged-ab3.par` (Adams-Bashforth 3, proposed schedule) and `charged-current.par` (`use_tapered_grids = yes`, `use_time_interpolation = no`, `use_final_copy = yes`).
