#include <assert.h>
#include <stdlib.h>

#include <iostream>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctki_GHExtensions.h"
#include "cctki_ScheduleBindings.h"
#include "cctki_WarnLevel.h"

#include "carpet.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Initialise.cc,v 1.45 2004/05/21 18:16:23 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_Carpet_Initialise_cc);
}



namespace Carpet {
  
  using namespace std;
  
  
  
  int Initialise (tFleshConfig* fc)
  {
    DECLARE_CCTK_PARAMETERS;
    
    // Initialise stuff
    const int convlev = 0;
    cGH* const cgh = CCTK_SetupGH (fc, convlev);
    CCTKi_AddGH (fc, convlev, cgh);
    
    // Delay checkpoint until MPI has been initialised
    Waypoint ("Starting initialisation");
    
    // Initialise stuff
    cgh->cctk_iteration = 0;
    global_time = cctk_initial_time;
    delta_time = 1.0;
    cgh->cctk_time = global_time;
    cgh->cctk_delta_time = delta_time;
    
    do_global_mode = true;
    do_meta_mode = true;
    
    // Enable storage and communtication
    CCTKi_ScheduleGHInit (cgh);
    
    // Initialise stuff
    CCTKi_InitGHExtensions (cgh);
    
    
    
    BEGIN_MGLEVEL_LOOP(cgh) {
      do_global_mode = true;
      do_meta_mode = mglevel==mglevels-1;
      
      // Register coordinates
      Checkpoint ("Scheduling CCTK_WRAGH");
      CCTK_ScheduleTraverse ("CCTK_WRAGH", cgh, CallFunction);
      
      // Check parameters
      Checkpoint ("Scheduling PARAMCHECK");
      CCTK_ScheduleTraverse ("CCTK_PARAMCHECK", cgh, CallFunction);
      CCTKi_FinaliseParamWarn();
    } END_MGLEVEL_LOOP;
    
    
    
    if (fc->recovered) {
      // if recovering
      
      
      
      for (int rl=0; rl<reflevels; ++rl) {
        BEGIN_MGLEVEL_LOOP(cgh) {
          enter_level_mode (cgh, rl);
          do_global_mode = reflevel==0;
          do_meta_mode = do_global_mode && mglevel==mglevels-1;
          
          cgh->cctk_time = global_time;
          
          Waypoint ("Recovering I at iteration %d time %g%s%s",
                    cgh->cctk_iteration, (double)cgh->cctk_time,
                    (do_global_mode ? " (global)" : ""),
                    (do_meta_mode ? " (meta)" : ""));
          
          // Set up the grids
          Checkpoint ("Scheduling BASEGRID");
          CCTK_ScheduleTraverse ("CCTK_BASEGRID", cgh, CallFunction);
          
          // Recover
          Checkpoint ("Scheduling RECOVER_VARIABLES");
          CCTK_ScheduleTraverse ("CCTK_RECOVER_VARIABLES", cgh, CallFunction);
          
          leave_level_mode (cgh);
        } END_MGLEVEL_LOOP;

        {
          const int ml=0;
          enter_global_mode (cgh, ml);
          enter_level_mode (cgh, rl);
          
          // Regrid
          Checkpoint ("Regrid");
          Regrid (cgh);
          
          // Postregrid
          Checkpoint ("Scheduling POSTREGRID");
          CCTK_ScheduleTraverse ("CCTK_POSTREGRID", cgh, CallFunction);
          
          leave_level_mode (cgh);
          leave_global_mode (cgh);
        } // ml
        
      } // for rl
      
      
      
      for (int rl=0; rl<reflevels; ++rl) {
        BEGIN_MGLEVEL_LOOP(cgh) {
          enter_level_mode (cgh, rl);
          do_global_mode = reflevel==reflevels-1;
          do_meta_mode = do_global_mode && mglevel==mglevels-1;
          
          Waypoint ("Recovering II at iteration %d time %g%s%s",
                    cgh->cctk_iteration, (double)cgh->cctk_time,
                    (do_global_mode ? " (global)" : ""),
                    (do_meta_mode ? " (meta)" : ""));
          
          // Post recover
          Checkpoint ("Scheduling POST_RECOVER_VARIABLES");
          CCTK_ScheduleTraverse
            ("CCTK_POST_RECOVER_VARIABLES", cgh, CallFunction);
          
          // Checking
          PoisonCheck (cgh, alltimes);
          CheckChecksums (cgh, allbutcurrenttime);
          
          leave_level_mode (cgh);
        } END_MGLEVEL_LOOP;
      } // for rl
      
      
      
    } else {
      // if not recovering
      
      
      
      for (int rl=0; rl<reflevels; ++rl) {
        BEGIN_MGLEVEL_LOOP(cgh) {
          enter_level_mode (cgh, rl);
          do_global_mode = reflevel==0;
          do_meta_mode = do_global_mode && mglevel==mglevels-1;
          
          cgh->cctk_time = global_time;
          
          Waypoint ("Initialisation I at iteration %d time %g%s%s",
                    cgh->cctk_iteration, (double)cgh->cctk_time,
                    (do_global_mode ? " (global)" : ""),
                    (do_meta_mode ? " (meta)" : ""));
          
          // Checking
          Poison (cgh, alltimes);
          
          // Set up the grids
          Checkpoint ("Scheduling BASEGRID");
          CCTK_ScheduleTraverse ("CCTK_BASEGRID", cgh, CallFunction);
          
          const int num_tl = init_each_timelevel ? 3 : 1;
          
          // Rewind
          for (int m=0; m<maps; ++m) {
            vtt.at(m)->set_delta
              (reflevel, mglevel, - vtt.at(m)->get_delta (reflevel, mglevel));
            FlipTimeLevels (cgh);
            for (int tl=0; tl<num_tl; ++tl) {
              vtt.at(m)->advance_time (reflevel, mglevel);
              CycleTimeLevels (cgh);
            }
            vtt.at(m)->set_delta
              (reflevel, mglevel, - vtt.at(m)->get_delta (reflevel, mglevel));
            FlipTimeLevels (cgh);
          }
          
          const bool outer_do_global_mode = do_global_mode;
          for (int tl=num_tl-1; tl>=0; --tl) {
            do_global_mode = outer_do_global_mode && tl==0;
            
            // Advance times
            for (int m=0; m<maps; ++m) {
              vtt.at(m)->advance_time (reflevel, mglevel);
            }
            cgh->cctk_time
              = global_time - tl * delta_time * mglevelfact / reflevelfact;
            CycleTimeLevels (cgh);
            
            // Set up the initial data
            Checkpoint ("Scheduling INITIAL");
            CCTK_ScheduleTraverse ("CCTK_INITIAL", cgh, CallFunction);
            
          } // for tl
          do_global_mode = outer_do_global_mode;
          
          // Checking
          PoisonCheck (cgh, currenttime);
          
          leave_level_mode (cgh);
        } END_MGLEVEL_LOOP;
        
        // Regrid
        {
          const int ml=0;
          enter_global_mode (cgh, ml);
          enter_level_mode (cgh, rl);
          
          // Regrid
          Checkpoint ("Regrid");
          Regrid (cgh);
          
          // Postregrid
          Checkpoint ("Scheduling POSTREGRID");
          CCTK_ScheduleTraverse ("CCTK_POSTREGRID", cgh, CallFunction);
          
          leave_level_mode (cgh);
          leave_global_mode (cgh);
        } // ml
        
      } // for rl
      
      
      
      for (int rl=reflevels-1; rl>=0; --rl) {
        BEGIN_MGLEVEL_LOOP(cgh) {
          enter_level_mode (cgh, rl);
          
          Waypoint ("Initialisation/Restrict at iteration %d time %g",
                    cgh->cctk_iteration, (double)cgh->cctk_time);
          
          // Restrict
          Restrict (cgh);
          
          leave_level_mode (cgh);
        } END_MGLEVEL_LOOP;
      } // for rl
      
      
      
      for (int rl=0; rl<reflevels; ++rl) {
        BEGIN_MGLEVEL_LOOP(cgh) {
          enter_level_mode (cgh, rl);
          do_global_mode = reflevel==reflevels-1;
          do_meta_mode = do_global_mode && mglevel==mglevels-1;
          
          Waypoint ("Initialisation II at iteration %d time %g%s%s",
                    cgh->cctk_iteration, (double)cgh->cctk_time,
                    (do_global_mode ? " (global)" : ""),
                    (do_meta_mode ? " (meta)" : ""));
          
          Checkpoint ("Scheduling POSTRESTRICTINITIAL");
          CCTK_ScheduleTraverse
            ("CCTK_POSTRESTRICTINITIAL", cgh, CallFunction);
          
          // Postinitial
          Checkpoint ("Scheduling POSTINITIAL");
          CCTK_ScheduleTraverse ("CCTK_POSTINITIAL", cgh, CallFunction);
          
          // Poststep
          Checkpoint ("Scheduling POSTSTEP");
          CCTK_ScheduleTraverse ("CCTK_POSTSTEP", cgh, CallFunction);
          
          // Checking
          PoisonCheck (cgh, alltimes);
          CheckChecksums (cgh, allbutcurrenttime);
          
          leave_level_mode (cgh);
        } END_MGLEVEL_LOOP;
      } // for rl
      
    
    
      if (init_3_timelevels) {
        // Use Scott Hawley's algorithm for getting two extra
        // timelevels of data
        Waypoint ("Initialising three timelevels");
      
        for (int rl=0; rl<reflevels; ++rl) {
          BEGIN_MGLEVEL_LOOP(cgh) {
            enter_level_mode (cgh, rl);
            do_global_mode = reflevel==0;
            do_meta_mode = do_global_mode && mglevel==mglevels-1;
          
            // Advance times
            for (int m=0; m<maps; ++m) {
              vtt.at(m)->advance_time (reflevel, mglevel);
            }
            cgh->cctk_time
              = global_time + delta_time * mglevelfact / reflevelfact;
            CycleTimeLevels (cgh);
          
            Waypoint ("Initialisation 3TL evolution I (a) (forwards) at iteration %d time %g%s%s",
                      cgh->cctk_iteration, (double)cgh->cctk_time,
                      (do_global_mode ? " (global)" : ""),
                      (do_meta_mode ? " (meta)" : ""));
          
            // Checking
            CalculateChecksums (cgh, allbutcurrenttime);
            Poison (cgh, currenttimebutnotifonly);
          
            // Evolve forward
            Checkpoint ("Scheduling PRESTEP");
            CCTK_ScheduleTraverse ("CCTK_PRESTEP", cgh, CallFunction);
            Checkpoint ("Scheduling EVOL");
            CCTK_ScheduleTraverse ("CCTK_EVOL", cgh, CallFunction);
          
            // Checking
            PoisonCheck (cgh, currenttime);
          
            leave_level_mode (cgh);
          } END_MGLEVEL_LOOP;
        } // for rl
      
        delta_time *= -1;
        for (int rl=0; rl<reflevels; ++rl) {
          BEGIN_MGLEVEL_LOOP(cgh) {
            enter_level_mode (cgh, rl);
            do_global_mode = reflevel==0;
            do_meta_mode = do_global_mode && mglevel==mglevels-1;
          
            // Flip time levels
            Waypoint ("Flipping timelevels");
            FlipTimeLevels (cgh);
          
            cgh->cctk_time
              = global_time + delta_time * mglevelfact / reflevelfact;
          
            leave_level_mode (cgh);
          } END_MGLEVEL_LOOP;
        } // for rl
      
        for (int rl=0; rl<reflevels; ++rl) {
          BEGIN_MGLEVEL_LOOP(cgh) {
            enter_level_mode (cgh, rl);
            do_global_mode = reflevel==0;
            do_meta_mode = do_global_mode && mglevel==mglevels-1;
          
            Waypoint ("Initialisation 3TL evolution I (b) (backwards) at iteration %d time %g%s%s",
                      cgh->cctk_iteration, (double)cgh->cctk_time,
                      (do_global_mode ? " (global)" : ""),
                      (do_meta_mode ? " (meta)" : ""));
          
            // Checking
            CalculateChecksums (cgh, allbutcurrenttime);
            Poison (cgh, currenttimebutnotifonly);
          
            // Evolve backward
            Checkpoint ("Scheduling PRESTEP");
            CCTK_ScheduleTraverse ("CCTK_PRESTEP", cgh, CallFunction);
            Checkpoint ("Scheduling EVOL");
            CCTK_ScheduleTraverse ("CCTK_EVOL", cgh, CallFunction);
          
            // Checking
            PoisonCheck (cgh, alltimes);
          
            leave_level_mode (cgh);
          } END_MGLEVEL_LOOP;
        } // for rl
      
        Waypoint ("Hourglass structure in place");
      
        // Evolve each level "backwards" one more timestep
        // Starting with the finest level and proceeding to the coarsest
        for (int rl=reflevels-1; rl>=0; --rl) {
          BEGIN_MGLEVEL_LOOP(cgh) {
            enter_level_mode (cgh, rl);
            do_global_mode = reflevel==0;
            do_meta_mode = do_global_mode && mglevel==mglevels-1;
            
            Waypoint ("Initialisation 3TL evolution II (b) (backwards) at iteration %d time %g%s%s",
                      cgh->cctk_iteration, (double)cgh->cctk_time,
                      (do_global_mode ? " (global)" : ""),
                      (do_meta_mode ? " (meta)" : ""));
            
            // Restrict
            Restrict (cgh);
            
            Checkpoint ("Scheduling POSTRESTRICT");
            CCTK_ScheduleTraverse ("CCTK_POSTRESTRICT", cgh, CallFunction);
            
            // Poststep
            Checkpoint ("Scheduling POSTSTEP");
            CCTK_ScheduleTraverse ("CCTK_POSTSTEP", cgh, CallFunction);
          
            // Checking
            PoisonCheck (cgh, alltimes);
          
            // Advance times
            for (int m=0; m<maps; ++m) {
              vtt.at(m)->advance_time (reflevel, mglevel);
            }
            cgh->cctk_time
              = global_time + 2 * delta_time * mglevelfact / reflevelfact;
            CycleTimeLevels (cgh);
          
            Waypoint ("Initialisation 3TL evolution I (c) (backwards) at iteration %d time %g%s%s",
                      cgh->cctk_iteration, (double)cgh->cctk_time,
                      (do_global_mode ? " (global)" : ""),
                      (do_meta_mode ? " (meta)" : ""));
          
            // Checking
            CalculateChecksums (cgh, allbutcurrenttime);
            Poison (cgh, currenttimebutnotifonly);
          
            // Evolve backward
            Checkpoint ("Scheduling PRESTEP");
            CCTK_ScheduleTraverse ("CCTK_PRESTEP", cgh, CallFunction);
            Checkpoint ("Scheduling EVOL");
            CCTK_ScheduleTraverse ("CCTK_EVOL", cgh, CallFunction);
            Checkpoint ("Scheduling POSTSTEP");
            CCTK_ScheduleTraverse ("CCTK_POSTSTEP", cgh, CallFunction);
          
            // Checking
            PoisonCheck (cgh, alltimes);
          
            leave_level_mode (cgh);
          } END_MGLEVEL_LOOP;
        } // for rl
      
        delta_time *= -1;
        for (int rl=0; rl<reflevels; ++rl) {
          BEGIN_MGLEVEL_LOOP(cgh) {
            enter_level_mode (cgh, rl);
            do_global_mode = reflevel==0;
            do_meta_mode = do_global_mode && mglevel==mglevels-1;
          
            // Flip time levels back
            Waypoint ("Flipping timelevels back");
            FlipTimeLevels (cgh);
          
            // Invert level times back
            for (int m=0; m<maps; ++m) {
              vtt.at(m)->set_delta
                (reflevel, mglevel, - vtt.at(m)->get_delta (reflevel, mglevel));
              vtt.at(m)->advance_time (reflevel, mglevel);
              vtt.at(m)->advance_time (reflevel, mglevel);
              vtt.at(m)->set_delta
                (reflevel, mglevel, - vtt.at(m)->get_delta (reflevel, mglevel));
            }
            cgh->cctk_time = global_time;
          
            leave_level_mode (cgh);
          } END_MGLEVEL_LOOP;
        } // for rl
      
        Waypoint ("Finished initialising three timelevels");
      
      } // if init_3_timelevels
    
    } // if not recovering
    
    
    
    for (int rl=0; rl<reflevels; ++rl) {
      BEGIN_MGLEVEL_LOOP(cgh) {
        enter_level_mode (cgh, rl);
        do_global_mode = reflevel==reflevels-1;
        do_meta_mode = do_global_mode && mglevel==mglevels-1;
        
        Waypoint ("Initialisation III at iteration %d time %g%s%s",
                  cgh->cctk_iteration, (double)cgh->cctk_time,
                  (do_global_mode ? " (global)" : ""),
                  (do_meta_mode ? " (meta)" : ""));
        
        // Checkpoint
        Checkpoint ("Scheduling CPINITIAL");
        CCTK_ScheduleTraverse ("CCTK_CPINITIAL", cgh, CallFunction);
        
        // Analysis
        Checkpoint ("Scheduling ANALYSIS");
        CCTK_ScheduleTraverse ("CCTK_ANALYSIS", cgh, CallFunction);
        
        // Output
        Checkpoint ("OutputGH");
        CCTK_OutputGH (cgh);
        
        // Checking
        PoisonCheck (cgh, alltimes);
        CheckChecksums (cgh, allbutcurrenttime);
        
        leave_level_mode (cgh);
      } END_MGLEVEL_LOOP;
    } // for rl
    
    
      
    Waypoint ("Done with initialisation");
    
    return 0;
  }
  
} // namespace Carpet
