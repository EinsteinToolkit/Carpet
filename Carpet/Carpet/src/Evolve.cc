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
  static const char* rcsid = "$Header:$";
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
    if (iteration % maxreflevelfact != 0) {
      
      // Terminate only after complete coarse grid steps
      // TODO: once checkpointing works correctly, change this to
      // "after complete time steps"
      return false;
      
    } else if (terminate_next || CCTK_TerminationReached(cgh)) {
      
      // Terminate if someone or something said so
      term = true;
      
    } else {
      
      const bool term_iter = iteration >= cctk_itlast;
      const bool term_time
        = (delta_time > 0
           ? time >= cctk_final_time - 1.0e-8 * cgh->cctk_delta_time
           : time <= cctk_final_time - 1.0e-8 * cgh->cctk_delta_time);
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
      } else if (CCTK_Equals(terminate, "immediately")) {
        term = true;
      } else {
        CCTK_WARN (0, "Unsupported termination condition");
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
  
  
  static void AdvanceTime( cGH* cgh, CCTK_REAL initial_time );
  static bool Regrid( cGH* cgh );
  static void PostRegrid( cGH* cgh );
  static void EvolutionI( cGH* cgh );
  static void Evolution_Restrict( cGH* cgh );
  static void EvolutionII( cGH* cgh );
  
  int Evolve (tFleshConfig* fc)
  {
    DECLARE_CCTK_PARAMETERS;
    
    Waypoint ("Starting evolution loop");
    
    const int convlev = 0;
    cGH* cgh = fc->GH[convlev];
    
    // Main loop
    while (! do_terminate(cgh, cgh->cctk_time, cgh->cctk_iteration)) {

      AdvanceTime( cgh, cctk_initial_time );

      if ((cgh->cctk_iteration-1)
          % (maxreflevelfact / ipow(reffact, reflevels-1)) == 0) {
        Waypoint ("Evolving iteration %d at t=%g",
                  cgh->cctk_iteration, (double)cgh->cctk_time);
      }
  
      if( Regrid( cgh ) )
        PostRegrid( cgh );
      
      EvolutionI( cgh );
      
      Evolution_Restrict( cgh );
      
      EvolutionII( cgh );
    }
    
    Waypoint ("Done with evolution loop");
    
    return 0;
  }

  void AdvanceTime( cGH* cgh, CCTK_REAL initial_time )
  {
      ++cgh->cctk_iteration;
      global_time = initial_time
        + cgh->cctk_iteration * delta_time / maxreflevelfact;
      cgh->cctk_time = global_time;
  }

  bool Regrid( cGH* cgh )
  {
    bool did_regrid = false;

    for (int rl=0; rl<reflevels; ++rl) {
      const int ml=0;
      const int do_every = maxreflevelfact / ipow(reffact, rl);
      if ((cgh->cctk_iteration-1) % do_every == 0) {
        enter_global_mode (cgh, ml);
        enter_level_mode (cgh, rl);
        
        Checkpoint ("Regrid");
        did_regrid |= Regrid (cgh, false, true);
         
        leave_level_mode (cgh);
        leave_global_mode (cgh);
      }
    }
    return did_regrid;
  }

  void PostRegrid( cGH* cgh )
  {
    for (int rl=0; rl<reflevels; ++rl) {
      for (int ml=mglevels-1; ml>=0; --ml) {
        enter_global_mode (cgh, ml);
        enter_level_mode (cgh, rl);
        
        do_global_mode = reflevel==0;
        do_meta_mode = do_global_mode && mglevel==mglevels-1;
        
        Waypoint ("Postregrid at iteration %d time %g%s%s",
                  cgh->cctk_iteration, (double)cgh->cctk_time,
                  (do_global_mode ? " (global)" : ""),
                  (do_meta_mode ? " (meta)" : ""));
        
        Checkpoint ("Scheduling POSTREGRID");
        CCTK_ScheduleTraverse ("CCTK_POSTREGRID", cgh, CallFunction);
        
        leave_level_mode (cgh);
        leave_global_mode (cgh);
      }
    }
  }

  void EvolutionI( cGH* cgh )
  {
    for (int ml=mglevels-1; ml>=0; --ml) {
      
      bool have_done_global_mode = false;
      bool have_done_anything = false;
      
      for (int rl=0; rl<reflevels; ++rl) {
        const int do_every
          = ipow(mgfact, ml) * (maxreflevelfact / ipow(reffact, rl));
        if ((cgh->cctk_iteration-1) % do_every == 0) {
          enter_global_mode (cgh, ml);
          enter_level_mode (cgh, rl);
          
          do_global_mode = ! have_done_global_mode;
          do_meta_mode = do_global_mode && mglevel==mglevels-1;
          assert (! (have_done_global_mode && do_global_mode));
          have_done_global_mode |= do_global_mode;
          have_done_anything = true;
          
          // Advance times
          for (int m=0; m<maps; ++m) {
            vtt.at(m)->advance_time (reflevel, mglevel);
          }
          cgh->cctk_time = (global_time
                            - delta_time / maxreflevelfact
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
          
          leave_level_mode (cgh);
          leave_global_mode (cgh);
        }
      }
      
      if (have_done_anything)
        assert (have_done_global_mode);
        
    }
  }

  void Evolution_Restrict( cGH* cgh )
  {
    for (int ml=mglevels-1; ml>=0; --ml) {
      for (int rl=reflevels-1; rl>=0; --rl) {
        const int do_every
          = ipow(mgfact, ml) * (maxreflevelfact / ipow(reffact, rl));
        if (cgh->cctk_iteration % do_every == 0) {
          enter_global_mode (cgh, ml);
          enter_level_mode (cgh, rl);
          
          Waypoint ("Evolution/Restrict at iteration %d time %g",
                    cgh->cctk_iteration, (double)cgh->cctk_time);
          
          Restrict (cgh);
          
          leave_level_mode (cgh);
          leave_global_mode (cgh);
        }
      }
    }
  }

  void EvolutionII( cGH* cgh )
  {
    for (int ml=mglevels-1; ml>=0; --ml) {
      
      bool have_done_global_mode = false;
      bool have_done_anything = false;
      
      for (int rl=0; rl<reflevels; ++rl) {
        const int do_every
          = ipow(mgfact, ml) * (maxreflevelfact / ipow(reffact, rl));
        if (cgh->cctk_iteration % do_every == 0) {
          enter_global_mode (cgh, ml);
          enter_level_mode (cgh, rl);
            
          int finest_active_reflevel = -1;
          {
            for (int rl_=0; rl_<reflevels; ++rl_) {
              const int do_every_
                = ipow(mgfact, ml) * (maxreflevelfact / ipow(reffact, rl_));
              if (cgh->cctk_iteration % do_every_ == 0) {
                finest_active_reflevel = rl_;
              }
            }
            assert (finest_active_reflevel >= 0);
          }
          do_global_mode = rl == finest_active_reflevel;
          do_meta_mode = do_global_mode && mglevel==mglevels-1;
          assert (! (have_done_global_mode && do_global_mode));
          have_done_global_mode |= do_global_mode;
          have_done_anything = true;
          
          Waypoint ("Evolution II at iteration %d time %g%s%s",
                    cgh->cctk_iteration, (double)cgh->cctk_time,
                    (do_global_mode ? " (global)" : ""),
                    (do_meta_mode ? " (meta)" : ""));
          
          Checkpoint ("Scheduling POSTRESTRICT");
          CCTK_ScheduleTraverse ("CCTK_POSTRESTRICT", cgh, CallFunction);
          
          // Poststep
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
          
          leave_level_mode (cgh);
          leave_global_mode (cgh);
        }
      }
 
      if (have_done_anything)
        assert (have_done_global_mode);
        
    }
  }
  
} // namespace Carpet
