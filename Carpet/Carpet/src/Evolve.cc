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
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Evolve.cc,v 1.30 2004/01/13 13:51:19 hawke Exp $";
  CCTK_FILEVERSION(Carpet_Carpet_Evolve_cc);
}



namespace Carpet {
  
  using namespace std;
  
  
  
  static double initial_time;
  
  
  
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
      double runtime;
#ifdef HAVE_TIME_GETTIMEOFDAY
      // get the current time
      struct timeval tv;
      gettimeofday (&tv, 0);
      runtime = (tv.tv_sec + tv.tv_usec / 1e6) - initial_time;
#else
      runtime = 0;
#endif
      const bool term_runtime = (max_runtime > 0
                                 && runtime >= 60.0 * max_runtime);
      
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
    
    Waypoint ("starting Evolve...");
    
    const int convlev = 0;
    cGH* cgh = fc->GH[convlev];
    
#ifdef HAVE_TIME_GETTIMEOFDAY
    // get the starting time
    struct timeval tv;
    gettimeofday (&tv, 0);
    initial_time = tv.tv_sec + tv.tv_usec / 1e6;
#else
    initial_time = 0;
#endif
    
    int next_global_mode_iter_loop1 = 0;
    int next_global_mode_iter_loop2 = 0;
    int next_global_mode_iter_loop3 = 0;
    
    // Main loop
    while (! do_terminate(cgh, refleveltimes[0], cgh->cctk_iteration)) {
      
      // Advance time
      ++cgh->cctk_iteration;
      
      Waypoint ("Evolving iteration %d...", cgh->cctk_iteration);
      
      BEGIN_REFLEVEL_LOOP(cgh) {
	const int do_every = maxreflevelfact/reflevelfact;
	if ((cgh->cctk_iteration-1) % do_every == 0) {
	  
	  BEGIN_MGLEVEL_LOOP(cgh) {
	    const int do_every = mglevelfact * (maxreflevelfact/reflevelfact);
	    if ((cgh->cctk_iteration-1) % do_every == 0) {
	      
              do_global_mode
                = cgh->cctk_iteration >= next_global_mode_iter_loop1;
              next_global_mode_iter_loop1 = cgh->cctk_iteration + 1;
              
	      // Advance level times
	      tt->advance_time (reflevel, mglevel);
              cgh->cctk_time = (cctk_initial_time
                                + (tt->time (0, reflevel, mglevel)
                                   * cgh->cctk_delta_time));
	      
	      Waypoint ("%*sCurrent time is %g, delta is %g%s", 2*reflevel, "",
			cgh->cctk_time,
                        cgh->cctk_delta_time / cgh->cctk_timefac,
                        do_global_mode ? "   (global time)" : "");
	      
	      // Cycle time levels
	      CycleTimeLevels (cgh);
	      
	      // Checking
	      CalculateChecksums (cgh, allbutcurrenttime);
	      Poison (cgh, currenttimebutnotifonly);
	      
	      // Evolve
	      Waypoint ("%*sScheduling PRESTEP", 2*reflevel, "");
	      CCTK_ScheduleTraverse ("CCTK_PRESTEP", cgh, CallFunction);
	      Waypoint ("%*sScheduling EVOL", 2*reflevel, "");
	      CCTK_ScheduleTraverse ("CCTK_EVOL", cgh, CallFunction);
	      Waypoint ("%*sScheduling POSTSTEP", 2*reflevel, "");
	      CCTK_ScheduleTraverse ("CCTK_POSTSTEP", cgh, CallFunction);
              
	      // Checking
	      PoisonCheck (cgh, currenttime);
	      
	    }
	  } END_MGLEVEL_LOOP;
          
	}
      } END_REFLEVEL_LOOP;
      
      
      
      BEGIN_REVERSE_REFLEVEL_LOOP(cgh) {
	const int do_every = maxreflevelfact/reflevelfact;
	if (cgh->cctk_iteration % do_every == 0) {
	  
	  BEGIN_MGLEVEL_LOOP(cgh) {
	    const int do_every = mglevelfact * (maxreflevelfact/reflevelfact);
	    if (cgh->cctk_iteration % do_every == 0) {
	      
              do_global_mode
                = cgh->cctk_iteration >= next_global_mode_iter_loop2;
              next_global_mode_iter_loop2 = cgh->cctk_iteration + 1;
              
	      // Restrict
	      Waypoint ("%*sCurrent time is %g, delta is %g%s", 2*reflevel, "",
			cgh->cctk_time,
                        cgh->cctk_delta_time / cgh->cctk_timefac,
                        do_global_mode ? "   (global time)" : "");
	      Restrict (cgh);
              
	      Waypoint ("%*sScheduling POSTRESTRICT", 2*reflevel, "");
	      CCTK_ScheduleTraverse ("CCTK_POSTRESTRICT", cgh, CallFunction);
	      
	      // Checking
	      CalculateChecksums (cgh, currenttime);
	      
	    }
	  } END_MGLEVEL_LOOP;
	  
	}
      } END_REVERSE_REFLEVEL_LOOP;
      
      
      
      BEGIN_REFLEVEL_LOOP(cgh) {
	const int do_every = maxreflevelfact/reflevelfact;
	if (cgh->cctk_iteration % do_every == 0) {
	  
	  // Regrid
	  Waypoint ("%*sRegrid", 2*reflevel, "");
	  Regrid (cgh, reflevel+1, true);
	  
	  BEGIN_MGLEVEL_LOOP(cgh) {
	    const int do_every = mglevelfact * (maxreflevelfact/reflevelfact);
	    if (cgh->cctk_iteration % do_every == 0) {
	      
              do_global_mode
                = cgh->cctk_iteration >= next_global_mode_iter_loop3;
              next_global_mode_iter_loop3 = cgh->cctk_iteration + 1;
	      
	      Waypoint ("%*sCurrent time is %g, delta is %g%s", 2*reflevel, "",
			cgh->cctk_time,
                        cgh->cctk_delta_time / cgh->cctk_timefac,
                        do_global_mode ? "   (global time)" : "");
              
	      // Checkpoint
	      Waypoint ("%*sScheduling CHECKPOINT", 2*reflevel, "");
	      CCTK_ScheduleTraverse ("CCTK_CHECKPOINT", cgh, CallFunction);
	      
	      // Analysis
	      Waypoint ("%*sScheduling ANALYSIS", 2*reflevel, "");
	      CCTK_ScheduleTraverse ("CCTK_ANALYSIS", cgh, CallFunction);
	      
	      // Output
	      Waypoint ("%*sOutputGH", 2*reflevel, "");
	      CCTK_OutputGH (cgh);
	      
	      // Checking
	      CheckChecksums (cgh, alltimes);
	      
	    }
	  } END_MGLEVEL_LOOP;
	  
	}
      } END_REFLEVEL_LOOP;
      
    } // main loop
    
    Waypoint ("done with Evolve.");
    
    return 0;
  }
  
} // namespace Carpet
