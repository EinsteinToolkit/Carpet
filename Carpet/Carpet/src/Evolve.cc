#include <assert.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Termination.h"

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

#include "dist.hh"
#include "th.hh"

#include "carpet.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Evolve.cc,v 1.34 2004/02/03 16:48:07 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_Carpet_Evolve_cc);
}



namespace Carpet {
  
  using namespace std;
  
  
  
  static bool do_terminate (const cGH *cgh,
                            const CCTK_REAL time, const int iteration)
  {
    DECLARE_CCTK_PARAMETERS;
    
    bool term;
    
    // Early shortcut
    if (terminate_next || CCTK_TerminationReached(cgh)) {
      
      term = true;
      
    } else {
      
      const bool term_iter = iteration >= cctk_itlast;
      const bool term_time = ((cctk_initial_time < cctk_final_time
                              ? time >= cctk_final_time
                              : time <= cctk_final_time)
                              && iteration % maxreflevelfact == 0);
#ifdef HAVE_TIME_GETTIMEOFDAY
      // get the current time
      struct timeval tv;
      gettimeofday (&tv, 0);
      const double thetime = tv.tv_sec + tv.tv_usec / 1e6;
      
      static bool firsttime = true;
      static double initial_runtime;
      if (firsttime) {
        firsttime = false;
        initial_runtime = thetime;
      }
      
      const double runtime = thetime - initial_runtime;
      const bool term_runtime = (max_runtime > 0
                                 && runtime >= 60.0 * max_runtime);
#else
      const bool term_runtime = false;
#endif
      
      if (CCTK_Equals(terminate, "never")) {
        term = false;
      } else if (CCTK_Equals(terminate, "iteration")) {
        term = term_iter;
      } else if (CCTK_Equals(terminate, "time")) {
        term = term_time;
      } else if (CCTK_Equals(terminate, "runtime")) {
        term = term_runtime;
      } else if (CCTK_Equals(terminate, "any")) {
        term = term_iter || term_time || term_runtime;
      } else if (CCTK_Equals(terminate, "all")) {
        term = term_iter && term_time && term_runtime;
      } else if (CCTK_Equals(terminate, "either")) {
        term = term_iter || term_time;
      } else if (CCTK_Equals(terminate, "both")) {
        term = term_iter && term_time;
      } else {
        assert (0);
      }
      
    }
    
    {
      int local, global;
      local = term;
      MPI_Allreduce (&local, &global, 1, MPI_INT, MPI_LOR, dist::comm);
      term = global;
    }
    
    return term;
  }
  
  
  
  int Evolve (tFleshConfig* fc)
  {
    DECLARE_CCTK_PARAMETERS;
    
    Waypoint ("Starting evolution loop");
    
    const int convlev = 0;
    cGH* cgh = fc->GH[convlev];
    
    // Main loop
    while (! do_terminate(cgh, cgh->cctk_time, cgh->cctk_iteration)) {
      
      // Regrid
      for (int rl=0; rl<reflevels; ++rl) {
        Checkpoint ("Regrid");
        Regrid (cgh, rl, rl+1, true);
      }
      
      // Advance time
      ++cgh->cctk_iteration;
      global_time += delta_time / maxreflevelfact;
      cgh->cctk_time = global_time;
      Waypoint ("Evolving iteration %d at t=%g",
                cgh->cctk_iteration, (double)cgh->cctk_time);
      
      bool have_done_global_mode = false;
      bool have_done_anything = false;
      for (int rl=0; rl<reflevels; ++rl) {
        BEGIN_MGLEVEL_LOOP(cgh) {
          enter_level_mode (cgh, rl);
          const int do_every = mglevelfact * (maxreflevelfact/reflevelfact);
          if ((cgh->cctk_iteration-1) % do_every == 0) {
            int coarsest_reflevel = maxreflevels;
            while (coarsest_reflevel > 0
                   && (cgh->cctk_iteration-1) % (mglevelfact * (maxreflevelfact / ipow(reffact, coarsest_reflevel-1))) == 0) {
              --coarsest_reflevel;
            }
            do_global_mode = reflevel==coarsest_reflevel;
            do_meta_mode = do_global_mode && mglevel==mglevels-1;
            if (do_global_mode) assert (! have_done_global_mode);
            have_done_global_mode = have_done_global_mode || do_global_mode;
            have_done_anything = true;
            
            // Advance times
            for (int m=0; m<maps; ++m) {
              vtt[m]->advance_time (reflevel, mglevel);
            }
            cgh->cctk_time = (global_time
                              - delta_time * mglevelfact / maxreflevelfact
                              + delta_time * mglevelfact / reflevelfact);
            CycleTimeLevels (cgh);
	    
            Waypoint ("Evolution I at iteration %d time %g%s%s",
                      cgh->cctk_iteration, (double)cgh->cctk_time,
                      (do_global_mode ? " (global)" : ""),
                      (do_meta_mode ? " (meta)" : ""));
            
            // Checking
            CalculateChecksums (cgh, allbutcurrenttime);
            Poison (cgh, currenttimebutnotifonly);
	    
            // Evolve
            Checkpoint ("Scheduling PRESTEP");
            CCTK_ScheduleTraverse ("CCTK_PRESTEP", cgh, CallFunction);
            Checkpoint ("Scheduling EVOL");
            CCTK_ScheduleTraverse ("CCTK_EVOL", cgh, CallFunction);
            
            // Checking
            PoisonCheck (cgh, currenttime);
            
          } // if do_every
          leave_level_mode (cgh);
        } END_MGLEVEL_LOOP;
      } // for rl
      if (have_done_anything) assert (have_done_global_mode);
      
      
      
      have_done_global_mode = false;
      have_done_anything = false;
      for (int rl=reflevels-1; rl>=0; --rl) {
        BEGIN_REVERSE_MGLEVEL_LOOP(cgh) {
          enter_level_mode (cgh, rl);
          const int do_every = mglevelfact * (maxreflevelfact/reflevelfact);
          if (cgh->cctk_iteration % do_every == 0) {
            int coarsest_reflevel = maxreflevels;
            while (coarsest_reflevel > 0
                   && cgh->cctk_iteration % (mglevelfact * (maxreflevelfact / ipow(reffact, coarsest_reflevel-1))) == 0) {
              --coarsest_reflevel;
            }
            do_global_mode = reflevel==coarsest_reflevel;
            do_meta_mode = do_global_mode && mglevel==mglevels-1;
            if (do_global_mode) assert (! have_done_global_mode);
            have_done_global_mode = have_done_global_mode || do_global_mode;
            have_done_anything = true;
            
            Waypoint ("Evolution II at iteration %d time %g%s%s",
                      cgh->cctk_iteration, (double)cgh->cctk_time,
                      (do_global_mode ? " (global)" : ""),
                      (do_meta_mode ? " (meta)" : ""));
            
            // Restrict
            Restrict (cgh);
            
            Checkpoint ("Scheduling POSTRESTRICT");
            CCTK_ScheduleTraverse ("CCTK_POSTRESTRICT", cgh, CallFunction);
            Checkpoint ("Scheduling POSTSTEP");
            CCTK_ScheduleTraverse ("CCTK_POSTSTEP", cgh, CallFunction);
	    
            // Checking
            PoisonCheck (cgh, currenttime);
            CalculateChecksums (cgh, currenttime);
            
            // Checkpoint
            Checkpoint ("Scheduling CHECKPOINT");
            CCTK_ScheduleTraverse ("CCTK_CHECKPOINT", cgh, CallFunction);
	    
            // Analysis
            Checkpoint ("Scheduling ANALYSIS");
            CCTK_ScheduleTraverse ("CCTK_ANALYSIS", cgh, CallFunction);
	    
            // Output
            Checkpoint ("OutputGH");
            CCTK_OutputGH (cgh);
	    
            // Checking
            CheckChecksums (cgh, alltimes);
	    
          } // if do_every
          leave_level_mode (cgh);
        } END_REVERSE_MGLEVEL_LOOP;
      } // for rl
      if (have_done_anything) assert (have_done_global_mode);
      
    } // main loop
    
    Waypoint ("Done with evolution loop");
    
    return 0;
  }
  
} // namespace Carpet
