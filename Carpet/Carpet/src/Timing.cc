#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

// IRIX wants this before <time.h>
#if HAVE_SYS_TYPES_H
#  include <sys/types.h>
#endif

#if TIME_WITH_SYS_TIME
#  include <sys/time.h>
#  include <time.h>
#else
#  if HAVE_SYS_TIME_H
#    include <sys/time.h>
#  elif HAVE_TIME_H
#    include <time.h>
#  endif
#endif

#if HAVE_UNISTD_H
#  include <unistd.h>
#endif

#include "defs.hh"
#include "dist.hh"

#include "carpet.hh"



namespace Carpet {
  
  using namespace std;
  
  
  
  // Return the current wall time
  static
  CCTK_REAL
  get_walltime ()
  {
#ifdef HAVE_TIME_GETTIMEOFDAY
    // get the current time
    struct timeval tv;
    gettimeofday (& tv, 0);
    return tv.tv_sec + tv.tv_usec / CCTK_REAL (1.0e6);
#else
    return CCTK_REAL (0.0);
#endif
  }
  
  
  
  // Calculate the number of updates for the current level
  static
  CCTK_REAL
  current_level_updates (cGH const * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    // Count the weighted number of grid points
    int num_grid_points = 0;
    for (int m = 0; m < maps; ++ m) {
      assert (reflevel >= 0);
      int const rl = reflevel;
      for (int c = 0; c < vhh.at(m)->components(rl); ++ c) {
        if (vhh.at(m)->is_local (rl, c)) {
          assert (mglevel >= 0);
          int const ml = mglevel;
          
          // Exterior of this region
          ibbox const ext = vdd.at(m)->boxes.at(ml).at(rl).at(c).exterior;
          // Outer boundaries
          b2vect const obs = xpose (vhh.at(m)->outer_boundary (rl, c));
          // Number of ghost zones
          i2vect const ghosts = vdd.at(m)->ghosts;
          // Computational domain: shrink the exterior by the number
          // of ghost zones if the face is not an outer boundary.
          // This ignores ghost zones, but takes buffer zones into
          // account.
          ibbox const domain =
            ext.expand (ivect (not obs[0]) * ghosts[0],
                        ivect (not obs[1]) * ghosts[1]);
          
          // Count the grid points
          num_grid_points += domain.size();
          
        } // if local
      }   // for c
    }     // for m
    
    // Take number of RHS evaluations per time step into account
    int const updates = num_grid_points * num_integrator_substeps;
    
    // Return result as CCTK_REAL to avoid potential overflows
    return updates;
  }
  
  
  
  // Time at which the evolution started
  CCTK_REAL initial_walltime;   // in seconds
  CCTK_REAL initial_phystime;
  
  // Counters for evolved grid points
  CCTK_REAL total_updates;
  
  
  
  // Initialise the timing variables (to be called during Initialise)
  void
  InitTimingVariables (cGH const * const cctkGH)
  {
    DECLARE_CCTK_ARGUMENTS;
    
    * grid_points_per_second = 0;
    * physical_time_per_hour = 0;
  }
  
  
  
  // Initialise the timing statistics (to be called at the beginning
  // of Evolve)
  void
  InitTiming (cGH const * const cctkGH)
  {
    initial_walltime = get_walltime();
    initial_phystime = cctkGH->cctk_time;
    
    total_updates = 0.0;
  }
  
  
  
  // Take a step on the current refinement and multigrid level (to be
  // called when EVOL is scheduled)
  void
  StepTiming (cGH const * const cctkGH)
  {
    total_updates += current_level_updates (cctkGH);
  }
  
  
  
  static
  void
  PrintUpdatesPerSecond (cGH const * const cctkGH)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    // Measure the elapsed time
    CCTK_REAL const elapsed_walltime = get_walltime() - initial_walltime;
    
    // Calculate updates per second
    CCTK_REAL const updates_per_second = total_updates / elapsed_walltime;
    
    * grid_points_per_second = updates_per_second;
    
    if (not silent) {
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Grid point updates per process per second: %g",
                  double (updates_per_second));
    }
    
#if 0
    CCTK_REAL const updates_per_second_2 = ipow (updates_per_second, 2);
    
    struct {
      CCTK_REAL ups, ups2;
    } local, global;
    local.ups  = updates_per_second;
    local.ups2 = updates_per_second_2;
    MPI_Allreduce (& local, & global, 2,
                   dist::datatype (global.ups), MPI_SUM, dist::comm());
    
    int const count = dist::size();
    CCTK_REAL const avg = global.ups / count;
    CCTK_REAL const stddev = sqrt (abs (global.ups2 - ipow (avg,2)) / count);
    
    CCTK_VInfo (CCTK_THORNSTRING,
                "Local updates per second:   %g", double (updates_per_second));
    CCTK_VInfo (CCTK_THORNSTRING,
                "Global updates per second:  %g", double (global.ups));
    
    if (verbose) {
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Average updates per second: %g", double (avg));
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Standard deviation:         %g", double (stddev));
    }
#endif
  }
  
  
  
  static
  void
  PrintPhysicalTimePerHour (cGH const * const cctkGH)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    // Measure the elapsed time
    CCTK_REAL const elapsed_walltime = get_walltime() - initial_walltime;
    
    // Calculate elapsed physical time
    CCTK_REAL const physical_time = cctkGH->cctk_time - initial_phystime;
    
    // Calculate physical time per hour
    * physical_time_per_hour =
      physical_time / elapsed_walltime * CCTK_REAL (3600.0);
    
    if (not silent) {
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Physical time per hour: %g",
                  double (* physical_time_per_hour));
    }
  }



  // Calculate and print some timing statistics (to be called at any
  // time)
  void
  PrintTimingStats (cGH const * const cctkGH)
  {
    PrintUpdatesPerSecond (cctkGH);
    PrintPhysicalTimePerHour (cctkGH);
  }
  
} // namespace Carpet
