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
  
  
  
  // Small number to avoid division by zero
  CCTK_REAL const eps = 1.0e-15;
  
  
  
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
                         CCTK_REAL & local_updates, CCTK_REAL & global_updates)
  {
    DECLARE_CCTK_PARAMETERS;
    
    // Count the weighted number of grid points
    // (int is not good enough for this calculation)
    CCTK_REAL local_num_grid_points = 0;
    CCTK_REAL global_num_grid_points = 0;
    for (int m = 0; m < maps; ++ m) {
      assert (reflevel >= 0);
      int const rl = reflevel;
      for (int c = 0; c < vhh.at(m)->components(rl); ++ c) {
        assert (mglevel >= 0);
        int const ml = mglevel;
        
        // Base region
        ibbox const ext = vhh.at(m)->extent(ml,rl,c);
        
        // Count the grid points
        CCTK_REAL const domainsize = ext.size();
        
        if (vhh.at(m)->is_local (rl, c)) {
          local_num_grid_points += domainsize;
        }
        global_num_grid_points += domainsize;
        
      } // for c
    }   // for m
    
    // Take number of RHS evaluations per time step into account
    static int int_steps = -1;
    if (int_steps == -1) {
      if (num_integrator_substeps != -1) {
        // if the user parameter is set, use it
        int_steps = num_integrator_substeps;
      } else if (CCTK_IsFunctionAliased ("MoLNumIntegratorSubsteps")) {
        // if there is an aliased function, use it
        int_steps = MoLNumIntegratorSubsteps ();
      } else {
        // use a sensible default, even if it is wrong -- it is better
        // to have timing information which is wrong by a constant
        // factor than to abort the code
        int_steps = 1;
      }
    }
    local_updates = local_num_grid_points * int_steps;
    global_updates = global_num_grid_points * int_steps;
  }
  
  
  
  // Time at which the evolution started
  CCTK_REAL initial_walltime;   // in seconds
  CCTK_REAL initial_phystime;
  
  
  
  // Last starting time
  enum timing_state_t { state_computing, state_communicating, state_io };
  timing_state_t timing_state = state_computing;
  CCTK_REAL time_start;
  
  
  
  // Initialise the timing variables (to be called before basegrid)
  void
  InitTimingStats (cGH const * const cctkGH)
  {
    DECLARE_CCTK_ARGUMENTS;
    
    * physical_time_per_hour = 0.0;
    
    * time_total         = 0.0;
    * time_computing     = 0.0;
    * time_communicating = 0.0;
    * time_io            = 0.0;
    
    * local_grid_points_per_second   = 0.0;
    * total_grid_points_per_second   = 0.0;
    * local_grid_point_updates_count = 0.0;
    * total_grid_point_updates_count = 0.0;
    
    * io_per_second              = 0.0;
    * io_bytes_per_second        = 0.0;
    * io_bytes_ascii_per_second  = 0.0;
    * io_bytes_binary_per_second = 0.0;
    * io_count                   = 0.0;
    * io_bytes_count             = 0.0;
    * io_bytes_ascii_count       = 0.0;
    * io_bytes_binary_count      = 0.0;
    
    * comm_per_second       = 0.0;
    * comm_bytes_per_second = 0.0;
    * comm_count            = 0.0;
    * comm_bytes_count      = 0.0;
    
    * grid_points_per_second   = 0.0;
    * grid_point_updates_count = 0.0;
  }
  
  
  
  // Begin timing (to be called after initialisation, just before the
  // main evolution begins)
  void
  BeginTiming (cGH const * const cctkGH)
  {
    DECLARE_CCTK_ARGUMENTS;
    
    initial_walltime = get_walltime();
    initial_phystime = cctkGH->cctk_time;
  }
  
  
  
  // Take a step on the current refinement and multigrid level (to be
  // called when EVOL is scheduled)
  void
  StepTiming (cGH const * const cctkGH)
  {
    DECLARE_CCTK_ARGUMENTS;
    
    assert (timing_state == state_computing);
    
    CCTK_REAL local_updates, global_updates;
    current_level_updates (cctkGH, local_updates, global_updates);
    
    * local_grid_point_updates_count += local_updates;
    * total_grid_point_updates_count += global_updates;
    
    * grid_point_updates_count = * local_grid_point_updates_count;
  }
  
  
  
  // Count some I/O (to be called from the I/O routine)
  void
  BeginTimingIO (cGH const * const cctkGH)
  {
    assert (timing_state == state_computing);
    timing_state = state_io;
    time_start = get_walltime();
  }
  
  void
  EndTimingIO (cGH const * const cctkGH,
               CCTK_REAL const files, CCTK_REAL const bytes,
               bool const is_binary)
  {
    DECLARE_CCTK_ARGUMENTS;
    
    assert (timing_state == state_io);
    timing_state = state_computing;
    CCTK_REAL const time_end = get_walltime();
    
    * time_io += time_end - time_start;
    
    * io_count += files;
    * io_bytes_count += bytes;
    * (is_binary ? io_bytes_binary_count : io_bytes_ascii_count) += bytes;
  }
  
  
  
  // Count some communication (to be called from the communication routine)
  void
  BeginTimingCommunication (cGH const * const cctkGH)
  {
    assert (timing_state == state_computing);
    timing_state = state_communicating;
    time_start = get_walltime();
  }
  
  void
  EndTimingCommunication (cGH const * const cctkGH,
                          CCTK_REAL const messages, CCTK_REAL const bytes)
  {
    DECLARE_CCTK_ARGUMENTS;
    
    assert (timing_state == state_communicating);
    timing_state = state_computing;
    CCTK_REAL const time_end = get_walltime();
    
    * time_communicating += time_end - time_start;
    
    * comm_count += messages;
    * comm_bytes_count += bytes;
  }
  
  
  
  static
  void
  UpdateTimes (cGH const * const cctkGH)
  {
    DECLARE_CCTK_ARGUMENTS;
    
    assert (timing_state == state_computing);
    
    // Measure the elapsed time
    * time_total = get_walltime() - initial_walltime;
    
    * time_computing = * time_total - (* time_communicating + * time_io);
  }
  
  
  
  static
  void
  UpdateUpdatesPerSecond (cGH const * const cctkGH)
  {
    DECLARE_CCTK_ARGUMENTS;
    
    // Calculate updates per second
    * local_grid_points_per_second =
      * local_grid_point_updates_count / max (* time_computing, eps);
    * total_grid_points_per_second =
      * total_grid_point_updates_count / max (* time_computing, eps);
    
    * grid_points_per_second = * local_grid_points_per_second;
  }
  
  
  
  static
  void
  UpdateIOStats (cGH const * const cctkGH)
  {
    DECLARE_CCTK_ARGUMENTS;
    
    * io_per_second =
      * io_count              / max (* time_io, eps);
    * io_bytes_per_second =
      * io_bytes_count        / max (* time_io, eps);
    * io_bytes_ascii_per_second =
      * io_bytes_ascii_count  / max (* time_io, eps);
    * io_bytes_binary_per_second =
      * io_bytes_binary_count / max (* time_io, eps);
  }
  
  
  
  static
  void
  UpdateCommunicationStats (cGH const * const cctkGH)
  {
    DECLARE_CCTK_ARGUMENTS;
    
    * comm_per_second =
      * comm_count       / max (* time_communicating, eps);
    * comm_bytes_per_second =
      * comm_bytes_count / max (* time_communicating, eps);
  }
  
  
  
  static
  void
  UpdatePhysicalTimePerHour (cGH const * const cctkGH)
  {
    DECLARE_CCTK_ARGUMENTS;
    
    // Calculate elapsed physical time
    CCTK_REAL const physical_time = cctkGH->cctk_time - initial_phystime;
    
    // Calculate physical time per hour
    * physical_time_per_hour = 3600.0 * physical_time / max (* time_total, eps);
  }
  
  
  
  // Calculate timing statistics (to be called before output)
  void
  UpdateTimingStats (cGH const * const cctkGH)
  {
    UpdateTimes (cctkGH);
    UpdateUpdatesPerSecond (cctkGH);
    UpdateIOStats (cctkGH);
    UpdateCommunicationStats (cctkGH);
    UpdatePhysicalTimePerHour (cctkGH);
  }
  
  
  
  static
  void
  PrintTimes (cGH const * const cctkGH)
  {
    DECLARE_CCTK_ARGUMENTS;
    
    CCTK_VInfo (CCTK_THORNSTRING,
                "Total run time: %g h", double (* time_total / 3600.0));
    CCTK_VInfo (CCTK_THORNSTRING,
                "(Comp, Comm, I/O) fractions = (%3.1f%%, %3.1f%%, %3.1f%%)",
                double (100.0 * * time_computing /
                        max (* time_total, eps)),
                double (100.0 * * time_communicating /
                        max (* time_total, eps)),
                double (100.0 * * time_io /
                        max (* time_total, eps)));
  }
  
  
  
  static
  void
  PrintUpdatesPerSecond (cGH const * const cctkGH)
  {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    CCTK_VInfo (CCTK_THORNSTRING,
                "This processor's grid point updates per second (local): %g",
                double (* local_grid_points_per_second));
    CCTK_VInfo (CCTK_THORNSTRING,
                "Overall grid point updates per second (total)         : %g",
                double (* total_grid_points_per_second));
    
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
                "Local updates per second:   %g",
                double (updates_per_second));
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
  PrintIOStats (cGH const * const cctkGH)
  {
    DECLARE_CCTK_ARGUMENTS;
    
    CCTK_VInfo (CCTK_THORNSTRING,
                "I/O operations per second:     %g",
                double (* io_per_second));
    CCTK_VInfo (CCTK_THORNSTRING,
                "I/O bytes per second:          %g",
                double (* io_bytes_per_second));
    CCTK_VInfo (CCTK_THORNSTRING,
                "I/O bytes per second (ASCII):  %g",
                double (* io_bytes_ascii_per_second));
    CCTK_VInfo (CCTK_THORNSTRING,
                "I/O bytes per second (binary): %g",
                double (* io_bytes_binary_per_second));
  }
  
  
  
  static
  void
  PrintCommunicationStats (cGH const * const cctkGH)
  {
    DECLARE_CCTK_ARGUMENTS;
    
    CCTK_VInfo (CCTK_THORNSTRING,
                "Communication operations per second: %g",
                double (* comm_per_second));
    CCTK_VInfo (CCTK_THORNSTRING,
                "Communicated bytes per second:       %g",
                double (* comm_bytes_per_second));
  }
  
  
  
  static
  void
  PrintPhysicalTimePerHour (cGH const * const cctkGH)
  {
    DECLARE_CCTK_ARGUMENTS;
    
    CCTK_VInfo (CCTK_THORNSTRING,
                "Physical time per hour: %g",
                double (* physical_time_per_hour));
  }
  
  
  
#define PRINTMEM(x) (memoryof(x) / 1.0e+6) << " MB"
  
  static
  void
  PrintMemoryStats (cGH const * const cctkGH)
  {
    DECLARE_CCTK_ARGUMENTS;
    
    cout << eol
         << "Memory statistics:" << eol
         << "   Grid hierarchy:" << eol;
    for (int m = 0; m < Carpet::maps; ++ m) {
      cout << "   gh[" << m << "]: " << PRINTMEM(*vhh.at(m)) << eol
           << "   dh[" << m << "]: " << PRINTMEM(*vdd.at(m)) << eol
           << "   th[" << m << "]: " << PRINTMEM(*vtt.at(m)) << eol;
    }
#if 0
    for (int g = 0; g < (int)arrdata.size(); ++ g) {
      if (CCTK_GroupTypeI(g) != CCTK_GF) {
        char * const groupname = CCTK_GroupName(g);
        for (int m = 0; m < (int)arrdata.at(g).size(); ++ m) {
          cout << "   Group " << groupname << ":" << eol
               << "   gh[" << m << "]: " << PRINTMEM(*arrdata.at(g).at(m).hh) << eol
               << "   dh[" << m << "]: " << PRINTMEM(*arrdata.at(g).at(m).dd) << eol
               << "   th[" << m << "]: " << PRINTMEM(*arrdata.at(g).at(m).tt) << eol;
          for (int v = 0; v < (int)arrdata.at(g).at(m).data.size(); ++ v) {
            char * const fullname = CCTK_FullName(CCTK_FirstVarIndexI(g)+v);
            cout << "   Variable " << fullname << ":" << eol
                 << "   ggf[" << m << "]: ";
            if (arrdata.at(g).at(m).data.at(v)) {
              cout << PRINTMEM(*arrdata.at(g).at(m).data.at(v));
            } else {
              cout << "<null>";
            }
            cout << eol;
            free (fullname);
          }
        }
        free (groupname);
      }
    }
#endif
    cout << endl;
  }



  // Calculate and print some timing statistics (to be called at any
  // time)
  void
  PrintTimingStats (cGH const * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (print_timestats_every > 0 and
        cctkGH->cctk_iteration % print_timestats_every == 0)
    {
      PrintTimes (cctkGH);
      PrintUpdatesPerSecond (cctkGH);
      PrintIOStats (cctkGH);
      PrintCommunicationStats (cctkGH);
      PrintPhysicalTimePerHour (cctkGH);
      PrintMemoryStats (cctkGH);
    }
  }
  
} // namespace Carpet
