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
  void
  current_level_updates (cGH const * const cctkGH,
                         int & local_updates, int & global_updates)
  {
    DECLARE_CCTK_PARAMETERS;
    
    // Count the weighted number of grid points
    int local_num_grid_points = 0;
    int global_num_grid_points = 0;
    for (int m = 0; m < maps; ++ m) {
      assert (reflevel >= 0);
      int const rl = reflevel;
      for (int c = 0; c < vhh.at(m)->components(rl); ++ c) {
        assert (mglevel >= 0);
        int const ml = mglevel;
        
        // Base region
        ibbox const ext = vhh.at(m)->extent(ml,rl,c);
        // Refinement boundaries
        b2vect const rbs = vhh.at(m)->refinement_boundaries (rl, c);
        // Number of buffer zones
        i2vect const buffers = vdd.at(m)->buffers;
        // Computational domain: Add the number of buffer zones to the
        // base extent.  This takes buffer zones into account and
        // ignores ghost zones.
        ibbox const domain =
          ext.expand (ivect (rbs[0]) * buffers[0], ivect (rbs[1]) * buffers[1]);
        
        // Count the grid points
        int const domainsize = domain.size();
        
        if (vhh.at(m)->is_local (rl, c)) {
          local_num_grid_points += domainsize;
        }
        global_num_grid_points += domainsize;
        
      } // for c
    }   // for m
    
    // Take number of RHS evaluations per time step into account
    local_updates = local_num_grid_points * num_integrator_substeps;
    global_updates = global_num_grid_points * num_integrator_substeps;
  }
  
  
  
  // Time at which the evolution started
  CCTK_REAL initial_walltime;   // in seconds
  CCTK_REAL initial_phystime;
  
  // Counters for evolved grid points
  CCTK_REAL total_local_updates;
  CCTK_REAL total_global_updates;
  
  
  
  // Initialise the timing variables (to be called during Initialise)
  void
  InitTimingVariables (cGH const * const cctkGH)
  {
    DECLARE_CCTK_ARGUMENTS;
    
    * physical_time_per_hour = 0.0;
    * local_grid_points_per_second = 0.0;
    * total_grid_points_per_second = 0.0;
    * grid_points_per_second = 0.0;
  }
  
  
  
  // Initialise the timing statistics (to be called at the beginning
  // of Evolve)
  void
  InitTiming (cGH const * const cctkGH)
  {
    initial_walltime = get_walltime();
    initial_phystime = cctkGH->cctk_time;
    
    total_local_updates = 0.0;
    total_global_updates = 0.0;
  }
  
  
  
  // Take a step on the current refinement and multigrid level (to be
  // called when EVOL is scheduled)
  void
  StepTiming (cGH const * const cctkGH)
  {
    int local_updates, global_updates;
    current_level_updates (cctkGH, local_updates, global_updates);
    total_local_updates = local_updates;
    total_global_updates = global_updates;
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
    CCTK_REAL const local_updates_per_second =
      total_local_updates / elapsed_walltime;
    CCTK_REAL const global_updates_per_second =
      total_global_updates / elapsed_walltime;
    
    * local_grid_points_per_second = local_updates_per_second;
    * total_grid_points_per_second = global_updates_per_second;
    *       grid_points_per_second = local_updates_per_second;
    
    if (not silent) {
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Local grid point updates per process per second: %g",
                  double (* local_grid_points_per_second));
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Total grid point updates per process per second: %g",
                  double (* total_grid_points_per_second));
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
