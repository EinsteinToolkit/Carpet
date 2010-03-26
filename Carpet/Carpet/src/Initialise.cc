#include <cassert>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <sstream>

#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctki_GHExtensions.h>
#include <cctki_ScheduleBindings.h>
#include <cctki_WarnLevel.h>

#include "carpet.hh"
#include "Timers.hh"



namespace Carpet {
  
  using namespace std;
  
  
  
  static void CallSetup (cGH * cctkGH);
  static void CallRecoverVariables (cGH * cctkGH);
#if 0
  static void CallRegridRecoverMeta (cGH * cctkGH);
#endif
  static void CallRegridRecoverLevel (cGH * cctkGH);
  static void CallRegridInitialMeta (cGH * cctkGH);
  static void CallRegridInitialLevel (cGH * cctkGH);
  static void CallPostRecoverVariables (cGH * cctkGH);
  static void CallInitial (cGH * cctkGH);
  static void CallRestrict (cGH * cctkGH);
  static void CallPostInitial (cGH * cctkGH);
  static void CallAnalysis (cGH * cctkGH);
  
  static void Initialise3tl (cGH * cctkGH);
  
  static void print_internal_data ();
  
  static void ScheduleTraverse
  (char const * where, char const * name, cGH * cctkGH);
  static void OutputGH (char const * where, cGH * cctkGH);
  
  
  
  int
  Initialise (tFleshConfig * const fc)
  {
    DECLARE_CCTK_PARAMETERS;
    
    int const convlev = 0;
    cGH * const cctkGH = CCTK_SetupGH (fc, convlev);
    CCTKi_AddGH (fc, convlev, cctkGH);
    
    do_global_mode = true;
    do_early_global_mode = true;
    do_late_global_mode = true;
    do_meta_mode = true;
    do_early_meta_mode = true;
    do_late_meta_mode = true;
    global_time = cctk_initial_time;
    delta_time = 1.0;
    for (int ml = 0; ml < mglevel; ++ ml) {
      assert (leveltimes.at(ml).size() == 1);
      leveltimes.at(ml).at(0) = global_time;
    }
    
    cctkGH->cctk_iteration = 0;
    cctkGH->cctk_time = global_time;
    cctkGH->cctk_delta_time = delta_time;
    
    static Timer timer ("Initialise");
    timer.start();
    
    // Delay checkpoint until MPI has been initialised
    Waypoint ("Starting initialisation");
    
    CCTKi_ScheduleGHInit (cctkGH); // Enable storage and communication
    GroupsStorageCheck (cctkGH);
    do_warn_about_storage = true;
    
    CCTKi_InitGHExtensions (cctkGH);
    
#if 0
    // Write grid structure to file
    for (int m=0; m<maps; ++m) {
      OutputGridStructure (cctkGH, m, vhh.at(m)->regions);
      OutputGridCoordinates (cctkGH, m, vhh.at(m)->regions);
    } // for m
#endif
    
    CallSetup (cctkGH);
    
    if (fc->recovered) {
      // Read data from a checkpoint file
      
      CallRecoverVariables (cctkGH);
      CallPostRecoverVariables (cctkGH);
      print_internal_data ();
      
    } else {
      // Calculate initial data
      
      CallInitial (cctkGH);
      CallRestrict (cctkGH);
      CallPostInitial (cctkGH);
      print_internal_data ();
      
      if (init_3_timelevels) {
        Initialise3tl (cctkGH);
      }
    }
    
    // Analyse initial data
    CallAnalysis (cctkGH);
    print_internal_data ();
    
    timer.stop();
    if (output_timers_every > 0) {
      TimerSet::writeData (cctkGH, timer_file);
    }
    
    Waypoint ("Done with initialisation");
    
    return 0;
  }
  
  
  
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  
  
  
  void
  CallSetup (cGH * const cctkGH)
  {
    char const * const where = "Initialise::CallSetup";
    static Timer timer (where);
    timer.start();
    
    BEGIN_MGLEVEL_LOOP(cctkGH) {
      do_early_global_mode = true;
      do_late_global_mode = true;
      do_global_mode = true;
      do_early_meta_mode = mglevel==mglevels-1;
      do_late_meta_mode = mglevel==0;
      do_meta_mode = do_early_meta_mode; // on first iteration, coarsest grid
      
      // Checking
      {
        int const rl=0;
        ENTER_LEVEL_MODE (cctkGH, rl) {
          Poison (cctkGH, alltimes, CCTK_ARRAY);
        } LEAVE_LEVEL_MODE;
      }
      
      // Register coordinates
      ScheduleTraverse (where, "CCTK_WRAGH", cctkGH);
      
      // Check parameters
      ScheduleTraverse (where, "CCTK_PARAMCHECK", cctkGH);
    } END_MGLEVEL_LOOP;
    
    Requirements::CheckRequirements (cctkGH);
    
    CCTKi_FinaliseParamWarn();
    
    timer.stop();
  }
  
  
  
  void
  CallRecoverVariables (cGH * const cctkGH)
  {
    char const * const where = "Initialise::CallRecoverVariables";
    static Timer timer (where);
    timer.start();
    
    DECLARE_CCTK_PARAMETERS;
    
    for (int rl=0; rl<reflevels; ++rl) {
      BEGIN_MGLEVEL_LOOP(cctkGH) {
        ENTER_LEVEL_MODE (cctkGH, rl) {
          do_early_global_mode = reflevel==0;
          do_late_global_mode = reflevel==reflevels-1;
          do_early_meta_mode = do_early_global_mode and mglevel==mglevels-1;
          do_late_meta_mode = do_late_global_mode and mglevel==0;
          do_global_mode = do_early_global_mode; // on first iteration, coarsest grid
          do_meta_mode = do_early_meta_mode; // on first iteration, coarsest grid
          
          // cctkGH->cctk_time = global_time;
          
          Waypoint ("Recover I at iteration %d time %g%s%s",
                    cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
                    (do_global_mode ? " (global)" : ""),
                    (do_meta_mode ? " (meta)" : ""));
          
          // Checking
          Poison (cctkGH, alltimes, CCTK_GF);
          
          // Set up the grids
          ScheduleTraverse (where, "CCTK_BASEGRID", cctkGH);
          
          // Recover
          ScheduleTraverse (where, "CCTK_RECOVER_VARIABLES", cctkGH);
          
          // Timing statistics
          // (do this here, after cctk_time has been recovered)
          if (do_early_global_mode) {
            InitTimingStats (cctkGH);
          }
          
          if (regrid_during_recovery) {
            CallRegridRecoverLevel (cctkGH);
          }
          
        } LEAVE_LEVEL_MODE;
      } END_MGLEVEL_LOOP;
    } // for rl
    
#if 0
    // TODO: Maybe do not checkpoint (nor read back) ghost and buffer
    // zones
    CallRegridRecoverMeta (cctkGH);
#endif
    
    timer.stop();
  }
  
  
  
  void
  CallPostRecoverVariables (cGH * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    char const * const where = "Initialise::CallPostRecoverVariables";
    static Timer timer (where);
    timer.start();
    
    for (int rl=0; rl<reflevels; ++rl) {
      BEGIN_MGLEVEL_LOOP(cctkGH) {
        ENTER_LEVEL_MODE (cctkGH, rl) {
          do_early_global_mode = reflevel==0;
          do_late_global_mode = reflevel==reflevels-1;
          do_early_meta_mode = do_early_global_mode and mglevel==mglevels-1;
          do_late_meta_mode = do_late_global_mode and mglevel==0;
          do_global_mode = do_early_global_mode; // on first iteration, coarsest grid
          do_meta_mode = do_early_meta_mode; // on first iteration, coarsest grid
          
          Waypoint ("Recovering II at iteration %d time %g%s%s",
                    cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
                    (do_global_mode ? " (global)" : ""),
                    (do_meta_mode ? " (meta)" : ""));
          
          int const num_tl = prolongation_order_time+1;
          
          bool const old_do_allow_past_timelevels = do_allow_past_timelevels;
          do_allow_past_timelevels = false;
          
          // Rewind times
          for (int m=0; m<maps; ++m) {
            CCTK_REAL const old_delta =
              vtt.AT(m)->get_delta (reflevel, mglevel);
            vtt.AT(m)->set_delta (reflevel, mglevel, - old_delta);
          }
          FlipTimeLevels (cctkGH);
          for (int tl=0; tl<num_tl; ++tl) {
            for (int m=0; m<maps; ++m) {
              vtt.AT(m)->advance_time (reflevel, mglevel);
            }
            CycleTimeLevels (cctkGH);
          }
          for (int m=0; m<maps; ++m) {
            CCTK_REAL const old_delta =
              vtt.AT(m)->get_delta (reflevel, mglevel);
            vtt.AT(m)->set_delta (reflevel, mglevel, - old_delta);
          }
          FlipTimeLevels (cctkGH);
          CCTK_REAL const old_cctk_time = cctkGH->cctk_time;
          cctkGH->cctk_time -=
            num_tl * (cctkGH->cctk_delta_time / cctkGH->cctk_timefac);
          
          for (int tl=0; tl<num_tl; ++tl) {
            
            // Advance times
            for (int m=0; m<maps; ++m) {
              vtt.AT(m)->advance_time (reflevel, mglevel);
            }
            CycleTimeLevels (cctkGH);
            cctkGH->cctk_time += cctkGH->cctk_delta_time / cctkGH->cctk_timefac;
            
            // Post recover variables
            ScheduleTraverse (where, "CCTK_POST_RECOVER_VARIABLES", cctkGH);
            
            // Checking
            PoisonCheck (cctkGH, currenttime);
            
          } // for tl
          cctkGH->cctk_time = old_cctk_time;
          
          do_allow_past_timelevels = old_do_allow_past_timelevels;
          
          CheckChecksums (cctkGH, allbutcurrenttime);
          
        } LEAVE_LEVEL_MODE;
      } END_MGLEVEL_LOOP;
    } // for rl
    
    timer.stop();
  }
  
  
  
  void
  CallInitial (cGH * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    char const * const where = "Initialise::CallInitial";
    static Timer timer (where);
    timer.start();
    
    if (not regrid_during_initialisation) {
      // Regrid once in the beginning
      CallRegridInitialMeta (cctkGH);
    }
    
    // Poison early, since grid functions may be initialised in global
    // loop-local mode, ane we must not overwrite them accidentally
    for (int rl=0; rl<reflevels; ++rl) {
      BEGIN_MGLEVEL_LOOP(cctkGH) {
        ENTER_LEVEL_MODE (cctkGH, rl) {
          
          // Checking
          Poison (cctkGH, alltimes, CCTK_GF);
          
        } LEAVE_LEVEL_MODE;
      } END_MGLEVEL_LOOP;
    } // for rl
    
    for (int rl=0; rl<reflevels; ++rl) {
      BEGIN_MGLEVEL_LOOP(cctkGH) {
        ENTER_LEVEL_MODE (cctkGH, rl) {
          do_early_global_mode = reflevel==0;
          do_late_global_mode = reflevel==reflevels-1;
          do_early_meta_mode = do_early_global_mode and mglevel==mglevels-1;
          do_late_meta_mode = do_late_global_mode and mglevel==0;
          do_global_mode = do_early_global_mode; // on first iteration, coarsest grid
          do_meta_mode = do_early_meta_mode; // on first iteration, coarsest grid
          
          // cctkGH->cctk_time = global_time;
          
          Waypoint ("Initialisation I at iteration %d time %g%s%s",
                    cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
                    (do_global_mode ? " (global)" : ""),
                    (do_meta_mode ? " (meta)" : ""));
          
          // Timing statistics
          if (do_early_global_mode) {
            InitTimingStats (cctkGH);
          }
          
          // Set up the grids
          ScheduleTraverse (where, "CCTK_BASEGRID", cctkGH);
          
          int const num_tl =
            init_each_timelevel ? prolongation_order_time+1 : 1;
          bool const old_do_allow_past_timelevels = do_allow_past_timelevels;
          do_allow_past_timelevels =
            not CCTK_EQUALS (initial_data_setup_method, "init_single_level");
          
          for (int m=0; m<maps; ++m) {
            vtt.at(m)->set_delta
              (reflevel, mglevel, - vtt.at(m)->get_delta (reflevel, mglevel));
            FlipTimeLevels (cctkGH);
            for (int tl=0; tl<num_tl; ++tl) {
              vtt.at(m)->advance_time (reflevel, mglevel);
              CycleTimeLevels (cctkGH);
            }
            vtt.at(m)->set_delta
              (reflevel, mglevel, - vtt.at(m)->get_delta (reflevel, mglevel));
            FlipTimeLevels (cctkGH);
          } // for m
          
          for (int tl=num_tl-1; tl>=0; --tl) {
            
            // Advance times
            for (int m=0; m<maps; ++m) {
              vtt.at(m)->advance_time (reflevel, mglevel);
            }
            cctkGH->cctk_time =
              (+ global_time
               - tl * delta_time * mglevelfact / timereflevelfact);
            CycleTimeLevels (cctkGH);
            
            // Set up the initial data
            ScheduleTraverse (where, "CCTK_INITIAL", cctkGH);
            ScheduleTraverse (where, "CCTK_POSTINITIAL", cctkGH);
            
            if (init_fill_timelevels) {
              assert (tl==0);
              FillTimeLevels (cctkGH);
            }
            
            // Checking
            PoisonCheck (cctkGH, currenttime);
            
          } // for tl
          
          do_allow_past_timelevels = old_do_allow_past_timelevels;
          
          if (regrid_during_initialisation and mglevel==0) {
            // Regrid after initialising each level
            CallRegridInitialLevel (cctkGH);
          }
          
        } LEAVE_LEVEL_MODE;
      } END_MGLEVEL_LOOP;
    } // for rl
    
    timer.stop();
  }
  
  
  
  void
  CallRestrict (cGH * const cctkGH)
  {
    for (int rl=reflevels-1; rl>=0; --rl) {
      BEGIN_MGLEVEL_LOOP(cctkGH) {
        ENTER_LEVEL_MODE (cctkGH, rl) {
          
          Waypoint ("Initialisation/Restrict at iteration %d time %g",
                    cctkGH->cctk_iteration, (double)cctkGH->cctk_time);
          
          Restrict (cctkGH);
          
        } LEAVE_LEVEL_MODE;
      } END_MGLEVEL_LOOP;
    } // for rl
  }
  
  
  
  void
  CallPostInitial (cGH * const cctkGH)
  {
    char const * const where = "Initialise::CallPostInitial";
    static Timer timer (where);
    timer.start();
    
    for (int rl=0; rl<reflevels; ++rl) {
      BEGIN_MGLEVEL_LOOP(cctkGH) {
        ENTER_LEVEL_MODE (cctkGH, rl) {
          do_early_global_mode = reflevel==0;
          do_late_global_mode = reflevel==reflevels-1;
          do_early_meta_mode = do_early_global_mode and mglevel==mglevels-1;
          do_late_meta_mode = do_late_global_mode and mglevel==0;
          do_global_mode = do_late_global_mode; // on last iteration, finest grid
          do_meta_mode = do_late_meta_mode; // on last iteration, finest grid
          
          Waypoint ("Initialisation II at iteration %d time %g%s%s",
                    cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
                    (do_global_mode ? " (global)" : ""),
                    (do_meta_mode ? " (meta)" : ""));
          
          if (reflevel < reflevels-1) {
            ScheduleTraverse (where, "CCTK_POSTRESTRICTINITIAL", cctkGH);
          }
          
          ScheduleTraverse (where, "CCTK_POSTPOSTINITIAL", cctkGH);
          ScheduleTraverse (where, "CCTK_POSTSTEP", cctkGH);
          
          PoisonCheck (cctkGH, alltimes);
          CheckChecksums (cctkGH, allbutcurrenttime);
          
        } LEAVE_LEVEL_MODE;
      } END_MGLEVEL_LOOP;
    } // for rl
    
    timer.stop();
  }
  
  
  
  void
  CallAnalysis (cGH * const cctkGH)
  {
    char const * const where = "Initialise::CallAnalysis";
    static Timer timer (where);
    timer.start();
    
    for (int rl=0; rl<reflevels; ++rl) {
      BEGIN_MGLEVEL_LOOP(cctkGH) {
        ENTER_LEVEL_MODE (cctkGH, rl) {
          do_early_global_mode = reflevel==0;
          do_late_global_mode = reflevel==reflevels-1;
          do_early_meta_mode = do_early_global_mode and mglevel==mglevels-1;
          do_late_meta_mode = do_late_global_mode and mglevel==0;
          do_global_mode = do_late_global_mode; // on last iteration, finest grid
          do_meta_mode = do_late_meta_mode; // on last iteration, finest grid
          
          Waypoint ("Initialisation III at iteration %d time %g%s%s",
                    cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
                    (do_global_mode ? " (global)" : ""),
                    (do_meta_mode ? " (meta)" : ""));
          
          int const do_every =
            ipow(mgfact, mglevel) * (maxtimereflevelfact / timereffacts.at(rl));
          if (cctkGH->cctk_iteration % do_every == 0)
          {
            // Checkpoint
            ScheduleTraverse (where, "CCTK_CPINITIAL", cctkGH);
            
            // Analysis
            ScheduleTraverse (where, "CCTK_ANALYSIS", cctkGH);
            
            if (do_late_global_mode) {
              // Timing statistics
              UpdateTimingStats (cctkGH);
            }
            
            // Output
            OutputGH (where, cctkGH);
            
            // Checking
            PoisonCheck (cctkGH, alltimes);
            CheckChecksums (cctkGH, allbutcurrenttime);
          } // if do_every
          
        } LEAVE_LEVEL_MODE;
      } END_MGLEVEL_LOOP;
    } // for rl
    
    timer.stop();
  }
  
  
  
#if 0
  void
  CallRegrid (cGH * const cctkGH,
              bool const callpreregrid,
              bool const callregrid,
              bool const callpostregrid,
              bool const regridinitial)
  {
    DECLARE_CCTK_PARAMETERS;
    
    assert (is_level_mode());
    
    bool const old_do_global_mode = do_global_mode;
    bool const old_do_meta_mode = do_meta_mode;
    do_global_mode = true;
    do_meta_mode = true;
    
    // Preregrid
    if (callpreregrid) {
      if (regridinitial) {
        Waypoint ("Preregridinitial at iteration %d time %g%s%s",
                  cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
                  (do_global_mode ? " (global)" : ""),
                  (do_meta_mode ? " (meta)" : ""));
        ScheduleTraverse (where, "CCTK_PREREGRIDINITIAL", cctkGH);
      } else {
        Waypoint ("Preregrid at iteration %d time %g%s%s",
                  cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
                  (do_global_mode ? " (global)" : ""),
                  (do_meta_mode ? " (meta)" : ""));
        ScheduleTraverse (where, "CCTK_PREREGRID", cctkGH);
      }
    }
    
    // Regrid
    bool did_regrid = false;
    if (callregrid) {
      Checkpoint ("Regrid");
      did_regrid = Regrid (cctkGH, true);
    }
    
    if (did_regrid or not callregrid) {
      BEGIN_META_MODE (cctkGH) {
        for (int rl=0; rl<reflevels; ++rl) {
          
          bool did_recompose = false;
          if (did_regrid) {
            did_recompose = Recompose (cctkGH, rl, prolongate_initial_data);
          }
          
          // Call postregridinitial only if initial data have already
          // been set up
          if (callpostregrid and (did_recompose or not callregrid)) {
            BEGIN_MGLEVEL_LOOP (cctkGH) {
              ENTER_LEVEL_MODE (cctkGH, rl) {
                do_global_mode = reflevel == reflevels - 1;
                do_meta_mode = do_global_mode and mglevel==mglevels-1;
                
                if (regridinitial) {
                  Waypoint ("Postregridinitial at iteration %d time %g%s%s",
                            cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
                            (do_global_mode ? " (global)" : ""),
                            (do_meta_mode ? " (meta)" : ""));
                } else {
                  Waypoint ("Postregrid at iteration %d time %g%s%s",
                            cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
                            (do_global_mode ? " (global)" : ""),
                            (do_meta_mode ? " (meta)" : ""));
                }
                
                int const num_tl =
                  regridinitial
                  ? (init_each_timelevel ? prolongation_order_time+1 : 1)
                  : prolongation_order_time+1;
                
                bool const old_do_allow_past_timelevels =
                  do_allow_past_timelevels;
                do_allow_past_timelevels = false;
                
                // Rewind times
                for (int m=0; m<maps; ++m) {
                  vtt.at(m)->set_delta
                    (reflevel, mglevel,
                     - vtt.at(m)->get_delta (reflevel, mglevel));
                  FlipTimeLevels (cctkGH);
                  for (int tl=0; tl<num_tl; ++tl) {
                    vtt.at(m)->advance_time (reflevel, mglevel);
                    CycleTimeLevels (cctkGH);
                  }
                  vtt.at(m)->set_delta
                    (reflevel, mglevel,
                     - vtt.at(m)->get_delta (reflevel, mglevel));
                  FlipTimeLevels (cctkGH);
                } // for m
                CCTK_REAL const old_cctk_time = cctkGH->cctk_time;
                cctkGH->cctk_time -=
                  num_tl * (cctkGH->cctk_delta_time / cctkGH->cctk_timefac);
                
                for (int tl=num_tl-1; tl>=0; --tl) {
                  
                  // Advance times
                  for (int m=0; m<maps; ++m) {
                    vtt.at(m)->advance_time (reflevel, mglevel);
                  }
                  CycleTimeLevels (cctkGH);
                  cctkGH->cctk_time +=
                    cctkGH->cctk_delta_time / cctkGH->cctk_timefac;
                  
                  // Postregrid
                  if (regridinitial) {
                    ScheduleTraverse (where, "CCTK_POSTREGRIDINITIAL", cctkGH);
                  } else {
                    ScheduleTraverse (where, "CCTK_POSTREGRID", cctkGH);
                  }
                  
                } // for tl
                cctkGH->cctk_time = old_cctk_time;
                
                do_allow_past_timelevels = old_do_allow_past_timelevels;
                
              } LEAVE_LEVEL_MODE;
            } END_MGLEVEL_LOOP;
          } // if did_recompose
          
        } // for rl
      } END_META_MODE;
    } // if did_regrid
    
    do_global_mode = old_do_global_mode;
    do_meta_mode = old_do_meta_mode;
  }
#endif
  
  
  
#if 0
  void
  CallRegridRecoverMeta (cGH * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    char const * const where = "Initialise::CallRegridRecoverMeta";
    static Timer timer (where);
    timer.start();
    
    assert (is_meta_mode());
    
    bool const old_do_global_mode = do_global_mode;
    bool const old_do_meta_mode = do_meta_mode;
    do_global_mode = true;
    do_meta_mode = true;
    
    for (int rl=0; rl<reflevels; ++rl) {
      
      BEGIN_MGLEVEL_LOOP (cctkGH) {
        ENTER_LEVEL_MODE (cctkGH, rl) {
          do_global_mode = reflevel == reflevels - 1;
          do_meta_mode = do_global_mode and mglevel==mglevels-1;
          
          Waypoint ("Postregrid at iteration %d time %g%s%s",
                    cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
                    (do_global_mode ? " (global)" : ""),
                    (do_meta_mode ? " (meta)" : ""));
          
          int const num_tl = prolongation_order_time + 1;
          
          bool const old_do_allow_past_timelevels = do_allow_past_timelevels;
          do_allow_past_timelevels = false;
          
          // Rewind times
          for (int m=0; m<maps; ++m) {
            vtt.at(m)->set_delta
              (reflevel, mglevel, - vtt.at(m)->get_delta (reflevel, mglevel));
            FlipTimeLevels (cctkGH);
            for (int tl=0; tl<num_tl; ++tl) {
              vtt.at(m)->advance_time (reflevel, mglevel);
              CycleTimeLevels (cctkGH);
            }
            vtt.at(m)->set_delta
              (reflevel, mglevel, - vtt.at(m)->get_delta (reflevel, mglevel));
            FlipTimeLevels (cctkGH);
          } // for m
          CCTK_REAL const old_cctk_time = cctkGH->cctk_time;
          cctkGH->cctk_time -=
            num_tl * (cctkGH->cctk_delta_time / cctkGH->cctk_timefac);
          
          for (int tl=num_tl-1; tl>=0; --tl) {
            
            // Advance times
            for (int m=0; m<maps; ++m) {
              vtt.at(m)->advance_time (reflevel, mglevel);
            }
            CycleTimeLevels (cctkGH);
            cctkGH->cctk_time +=
              cctkGH->cctk_delta_time / cctkGH->cctk_timefac;
            
            // Postregrid
            ScheduleTraverse (where, "CCTK_POSTREGRID", cctkGH);
              
          } // for tl
          cctkGH->cctk_time = old_cctk_time;
          
          do_allow_past_timelevels = old_do_allow_past_timelevels;
          
        } LEAVE_LEVEL_MODE;
      } END_MGLEVEL_LOOP;
    } // for rl
    
    do_global_mode = old_do_global_mode;
    do_meta_mode = old_do_meta_mode;
    
    timer.stop();
  }
#endif
  
  
  
  void
  CallRegridRecoverLevel (cGH * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    char const * const where = "Initialise::CallRegridRecoverLevel";
    static Timer timer (where);
    timer.start();
    
    CCTK_WARN (CCTK_WARN_ALERT,
               "Regridding in level mode after recovering is discouraged");
    
    assert (is_level_mode());
    
    bool const old_do_global_mode = do_global_mode;
    bool const old_do_early_global_mode = do_early_global_mode;
    bool const old_do_late_global_mode = do_late_global_mode;
    bool const old_do_meta_mode = do_meta_mode;
    bool const old_do_early_meta_mode = do_early_meta_mode;
    bool const old_do_late_meta_mode = do_late_meta_mode;
    do_global_mode = true;
    do_early_global_mode = true;
    do_late_global_mode = true;
    do_meta_mode = true;
    do_early_meta_mode = true;
    do_late_meta_mode = true;
    
    // Preregrid
    Waypoint ("Preregrid at iteration %d time %g%s%s",
              cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
              (do_global_mode ? " (global)" : ""),
              (do_meta_mode ? " (meta)" : ""));
    ScheduleTraverse (where, "CCTK_PREREGRID", cctkGH);
    
    // Regrid
    Checkpoint ("Regrid");
    bool const did_regrid = Regrid (cctkGH, true);
    
    if (did_regrid) {
      BEGIN_META_MODE (cctkGH) {
        for (int rl=0; rl<reflevels; ++rl) {
          
          bool did_recompose = false;
          if (did_regrid) {
            did_recompose = Recompose (cctkGH, rl, prolongate_initial_data);
          }
          
          if (did_recompose) {
            BEGIN_MGLEVEL_LOOP (cctkGH) {
              ENTER_LEVEL_MODE (cctkGH, rl) {
                do_early_global_mode = reflevel==0;
                do_late_global_mode = reflevel==reflevels-1;
                do_early_meta_mode = do_early_global_mode and mglevel==mglevels-1;
                do_late_meta_mode = do_late_global_mode and mglevel==0;
                do_global_mode = do_late_global_mode;
                do_meta_mode = do_late_meta_mode;
                
                Waypoint ("Postregrid at iteration %d time %g%s%s",
                          cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
                          (do_global_mode ? " (global)" : ""),
                          (do_meta_mode ? " (meta)" : ""));
                
                int const num_tl = prolongation_order_time + 1;
                
                bool const old_do_allow_past_timelevels =
                  do_allow_past_timelevels;
                do_allow_past_timelevels = false;
                
                // Rewind times
                for (int m=0; m<maps; ++m) {
                  vtt.at(m)->set_delta
                    (reflevel, mglevel,
                     - vtt.at(m)->get_delta (reflevel, mglevel));
                  FlipTimeLevels (cctkGH);
                  for (int tl=0; tl<num_tl; ++tl) {
                    vtt.at(m)->advance_time (reflevel, mglevel);
                    CycleTimeLevels (cctkGH);
                  }
                  vtt.at(m)->set_delta
                    (reflevel, mglevel,
                     - vtt.at(m)->get_delta (reflevel, mglevel));
                  FlipTimeLevels (cctkGH);
                } // for m
                CCTK_REAL const old_cctk_time = cctkGH->cctk_time;
                cctkGH->cctk_time -=
                  num_tl * (cctkGH->cctk_delta_time / cctkGH->cctk_timefac);
                
                for (int tl=num_tl-1; tl>=0; --tl) {
                  
                  // Advance times
                  for (int m=0; m<maps; ++m) {
                    vtt.at(m)->advance_time (reflevel, mglevel);
                  }
                  CycleTimeLevels (cctkGH);
                  cctkGH->cctk_time +=
                    cctkGH->cctk_delta_time / cctkGH->cctk_timefac;
                  
                  // Postregrid
                  ScheduleTraverse (where, "CCTK_POSTREGRID", cctkGH);
                  
                } // for tl
                cctkGH->cctk_time = old_cctk_time;
                
                do_allow_past_timelevels = old_do_allow_past_timelevels;
                
              } LEAVE_LEVEL_MODE;
            } END_MGLEVEL_LOOP;
          } // if did_recompose
          
        } // for rl
      } END_META_MODE;
    } // if did_regrid
    
    do_global_mode = old_do_global_mode;
    do_early_global_mode = old_do_early_global_mode;
    do_late_global_mode = old_do_late_global_mode;
    do_meta_mode = old_do_meta_mode;
    do_early_meta_mode = old_do_early_meta_mode;
    do_late_meta_mode = old_do_late_meta_mode;
    
    timer.stop();
  }
  
  
  
  void
  CallRegridInitialMeta (cGH * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    char const * const where = "Initialise::CallRegridInitialMeta";
    static Timer timer (where);
    timer.start();
    
    assert (is_meta_mode());
    
    bool const old_do_global_mode = do_global_mode;
    bool const old_do_early_global_mode = do_early_global_mode;
    bool const old_do_late_global_mode = do_late_global_mode;
    bool const old_do_meta_mode = do_meta_mode;
    bool const old_do_early_meta_mode = do_early_meta_mode;
    bool const old_do_late_meta_mode = do_late_meta_mode;
    do_global_mode = true;
    do_early_global_mode = true;
    do_late_global_mode = true;
    do_meta_mode = true;
    do_early_meta_mode = true;
    do_late_meta_mode = true;
    
    ENTER_GLOBAL_MODE (cctkGH, 0) {
      ENTER_LEVEL_MODE (cctkGH, 0) {
        
        // Preregrid
        Waypoint ("Preregridinitial at iteration %d time %g%s%s",
                  cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
                  (do_global_mode ? " (global)" : ""),
                  (do_meta_mode ? " (meta)" : ""));
        ScheduleTraverse (where, "CCTK_PREREGRIDINITIAL", cctkGH);
        
        // Regrid
        Checkpoint ("Regrid");
        bool const did_regrid = Regrid (cctkGH, true);
        
        if (did_regrid) {
          for (int rl=0; rl<reflevels; ++rl) {
            Recompose (cctkGH, rl, prolongate_initial_data);
          } // for rl
        } // if did_regrid
        
      } LEAVE_LEVEL_MODE;
    } LEAVE_GLOBAL_MODE;
    
    do_global_mode = old_do_global_mode;
    do_early_global_mode = old_do_early_global_mode;
    do_late_global_mode = old_do_late_global_mode;
    do_meta_mode = old_do_meta_mode;
    do_early_meta_mode = old_do_early_meta_mode;
    do_late_meta_mode = old_do_late_meta_mode;
    
    timer.stop();
  }
  
  
  
  void
  CallRegridInitialLevel (cGH * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    char const * const where = "Initialise::CallRegridInitialLevel";
    static Timer timer (where);
    timer.start();
    
    CCTK_WARN (CCTK_WARN_ALERT,
               "Regridding in level mode while initialising is discouraged");
    
    assert (is_level_mode());
    
    bool const old_do_global_mode = do_global_mode;
    bool const old_do_early_global_mode = do_early_global_mode;
    bool const old_do_late_global_mode = do_late_global_mode;
    bool const old_do_meta_mode = do_meta_mode;
    bool const old_do_early_meta_mode = do_early_meta_mode;
    bool const old_do_late_meta_mode = do_late_meta_mode;
    do_global_mode = true;
    do_early_global_mode = true;
    do_late_global_mode = true;
    do_meta_mode = true;
    do_early_meta_mode = true;
    do_late_meta_mode = true;
    
    // Preregrid
    Waypoint ("Preregridinitial at iteration %d time %g%s%s",
              cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
              (do_global_mode ? " (global)" : ""),
              (do_meta_mode ? " (meta)" : ""));
    ScheduleTraverse (where, "CCTK_PREREGRIDINITIAL", cctkGH);
    
    // Regrid
    Checkpoint ("Regrid");
    bool const did_regrid = Regrid (cctkGH, true);
    
    if (did_regrid) {
      BEGIN_META_MODE (cctkGH) {
        for (int rl=0; rl<reflevels; ++rl) {
          
          bool did_recompose = false;
          if (did_regrid) {
            did_recompose = Recompose (cctkGH, rl, prolongate_initial_data);
          }
          
          if (did_recompose) {
            BEGIN_MGLEVEL_LOOP (cctkGH) {
              ENTER_LEVEL_MODE (cctkGH, rl) {
                do_early_global_mode = reflevel==0;
                do_late_global_mode = reflevel==reflevels-1;
                do_early_meta_mode = do_early_global_mode and mglevel==mglevels-1;
                do_late_meta_mode = do_late_global_mode and mglevel==0;
                do_global_mode = do_late_global_mode;
                do_meta_mode = do_late_meta_mode;
                
                Waypoint ("Postregridinitial at iteration %d time %g%s%s",
                          cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
                          (do_global_mode ? " (global)" : ""),
                          (do_meta_mode ? " (meta)" : ""));
                
                int const num_tl =
                  init_each_timelevel ? prolongation_order_time+1 : 1;
                
                bool const old_do_allow_past_timelevels =
                  do_allow_past_timelevels;
                do_allow_past_timelevels = false;
                
                // Rewind times
                for (int m=0; m<maps; ++m) {
                  vtt.at(m)->set_delta
                    (reflevel, mglevel,
                     - vtt.at(m)->get_delta (reflevel, mglevel));
                  FlipTimeLevels (cctkGH);
                  for (int tl=0; tl<num_tl; ++tl) {
                    vtt.at(m)->advance_time (reflevel, mglevel);
                    CycleTimeLevels (cctkGH);
                  }
                  vtt.at(m)->set_delta
                    (reflevel, mglevel,
                     - vtt.at(m)->get_delta (reflevel, mglevel));
                  FlipTimeLevels (cctkGH);
                } // for m
                CCTK_REAL const old_cctk_time = cctkGH->cctk_time;
                cctkGH->cctk_time -=
                  num_tl * (cctkGH->cctk_delta_time / cctkGH->cctk_timefac);
                
                for (int tl=num_tl-1; tl>=0; --tl) {
                  
                  // Advance times
                  for (int m=0; m<maps; ++m) {
                    vtt.at(m)->advance_time (reflevel, mglevel);
                  }
                  CycleTimeLevels (cctkGH);
                  cctkGH->cctk_time +=
                    cctkGH->cctk_delta_time / cctkGH->cctk_timefac;
                  
                  // Postregrid
                  ScheduleTraverse (where, "CCTK_POSTREGRIDINITIAL", cctkGH);
                  
                } // for tl
                cctkGH->cctk_time = old_cctk_time;
                
                do_allow_past_timelevels = old_do_allow_past_timelevels;
                
              } LEAVE_LEVEL_MODE;
            } END_MGLEVEL_LOOP;
          } // if did_recompose
          
        } // for rl
      } END_META_MODE;
    } // if did_regrid
    
    do_global_mode = old_do_global_mode;
    do_early_global_mode = old_do_early_global_mode;
    do_late_global_mode = old_do_late_global_mode;
    do_meta_mode = old_do_meta_mode;
    do_early_meta_mode = old_do_early_meta_mode;
    do_late_meta_mode = old_do_late_meta_mode;
    
    timer.stop();
  }
  
  
  
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  
  
  
  // Use Scott Hawley's algorithm to get two extra timelevels of data
  
  static void initialise_3tl_advance_time (cGH * const cctkGH);
  static void initialise_3tl_evolve_Ia (cGH * const cctkGH);
  static void initialise_3tl_flip_timelevels (cGH * const cctkGH);
  static void initialise_3tl_evolve_Ib (cGH * const cctkGH);
  static void initialise_3tl_evolve_IIb (cGH * const cctkGH);
  static void initialise_3tl_advance_time_2 (cGH * const cctkGH);
  static void initialise_3tl_evolve_Ic (cGH * const cctkGH);
  static void initialise_3tl_reset_time (cGH * const cctkGH);
  
  void
  Initialise3tl (cGH * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    Waypoint ("Initialising three timelevels");
    
    // TODO: ensure that there are 3 timelevels
    assert (prolongation_order_time == 2);
    
    BEGIN_MGLEVEL_LOOP(cctkGH) {
      BEGIN_REFLEVEL_LOOP(cctkGH) {
        do_early_global_mode = reflevel==0;
        do_late_global_mode = reflevel==reflevels-1;
        do_early_meta_mode = do_early_global_mode and mglevel==mglevels-1;
        do_late_meta_mode = do_late_global_mode and mglevel==0;
        do_global_mode = do_early_global_mode;
        do_meta_mode = do_early_meta_mode;
        
        initialise_3tl_advance_time (cctkGH);
        initialise_3tl_evolve_Ia (cctkGH);
        initialise_3tl_flip_timelevels (cctkGH);
        initialise_3tl_evolve_Ib (cctkGH);
        initialise_3tl_flip_timelevels (cctkGH);
        
      } END_REFLEVEL_LOOP;
    } END_MGLEVEL_LOOP;
    
    Waypoint ("Hourglass structure in place");
    
    initialise_3tl_flip_timelevels (cctkGH);
    
    BEGIN_MGLEVEL_LOOP(cctkGH) {
      BEGIN_REVERSE_REFLEVEL_LOOP(cctkGH) {
        do_early_global_mode = reflevel==reflevels-1;
        do_late_global_mode = reflevel==0;
        do_early_meta_mode = do_early_global_mode and mglevel==mglevels-1;
        do_late_meta_mode = do_late_global_mode and mglevel==0;
        do_global_mode = do_early_global_mode;
        do_meta_mode = do_early_meta_mode;
        
        initialise_3tl_evolve_IIb (cctkGH);
        initialise_3tl_advance_time_2 (cctkGH);
        initialise_3tl_evolve_Ic (cctkGH);
        
      } END_REVERSE_REFLEVEL_LOOP;
    } END_MGLEVEL_LOOP;
    
    initialise_3tl_flip_timelevels (cctkGH);
    
    BEGIN_MGLEVEL_LOOP(cctkGH) {
      BEGIN_REVERSE_REFLEVEL_LOOP(cctkGH) {
        do_early_global_mode = reflevel==reflevels-1;
        do_late_global_mode = reflevel==0;
        do_early_meta_mode = do_early_global_mode and mglevel==mglevels-1;
        do_late_meta_mode = do_late_global_mode and mglevel==0;
        do_global_mode = do_early_global_mode;
        do_meta_mode = do_early_meta_mode;
        
        initialise_3tl_reset_time (cctkGH);
        
      } END_REVERSE_REFLEVEL_LOOP;
    } END_MGLEVEL_LOOP;
    
    Waypoint ("Finished initialising three timelevels");
  }
  
  
  
  void
  initialise_3tl_advance_time (cGH * const cctkGH)
  {
    Waypoint ("Advancing time");
    
    cctkGH->cctk_time
      = global_time + delta_time * mglevelfact / timereflevelfact;
    for (int m=0; m<maps; ++m) {
      vtt.at(m)->advance_time (reflevel, mglevel);
    }
    
    CycleTimeLevels (cctkGH);
  }
  
  void
  initialise_3tl_evolve_Ia (cGH * const cctkGH)
  {
    char const * const where = "Initialise3TL::EvolveIa";
    static Timer timer (where);
    timer.start();
    
    Waypoint ("Initialisation 3TL evolution I (a) (forwards) at iteration"
              " %d time %g%s%s",
              cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
              (do_global_mode ? " (global)" : ""),
              (do_meta_mode ? " (meta)" : ""));
    
    CalculateChecksums (cctkGH, allbutcurrenttime);
    Poison (cctkGH, currenttimebutnotifonly);
    
    // Evolve forward
    ScheduleTraverse (where, "CCTK_PRESTEP", cctkGH);
    ScheduleTraverse (where, "CCTK_EVOL", cctkGH);
    
    PoisonCheck (cctkGH, currenttime);
    
    timer.stop();
  }
  
  void
  initialise_3tl_flip_timelevels (cGH * const cctkGH)
  {
    Waypoint ("Flipping timelevels");
    
    BEGIN_META_MODE(cctkGH) {
      
      delta_time *= -1;
      
      BEGIN_MGLEVEL_LOOP(cctkGH) {
        BEGIN_REFLEVEL_LOOP (cctkGH) {
          
          cctkGH->cctk_time
            = global_time + delta_time * mglevelfact / timereflevelfact;
          
          FlipTimeLevels (cctkGH);
          
        } END_REFLEVEL_LOOP;
      } END_MGLEVEL_LOOP;
    } END_META_MODE;
  }
  
  void
  initialise_3tl_evolve_Ib (cGH * const cctkGH)
  {
    char const * const where = "Initialise3TL::EvolveIb";
    static Timer timer (where);
    timer.start();
    
    Waypoint ("Initialisation 3TL evolution I (b) (backwards) at iteration"
              " %d time %g%s%s",
              cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
              (do_global_mode ? " (global)" : ""),
              (do_meta_mode ? " (meta)" : ""));
    
    // Checking
    CalculateChecksums (cctkGH, allbutcurrenttime);
    Poison (cctkGH, currenttimebutnotifonly);
    
    // Evolve backward
    ScheduleTraverse (where, "CCTK_PRESTEP", cctkGH);
    ScheduleTraverse (where, "CCTK_EVOL", cctkGH);
    
    // Checking
    PoisonCheck (cctkGH, alltimes);
    
    timer.stop();
  }
  
  // Evolve backwards one more timestep
  // Starting with the finest level and proceeding to the coarsest
  void
  initialise_3tl_evolve_IIb (cGH * const cctkGH)
  {
    char const * const where = "Initialise3TL::EvolveIIb";
    static Timer timer (where);
    timer.start();
    
    Waypoint ("Initialisation 3TL evolution II (b) (backwards) at iteration"
              " %d time %g%s%s",
              cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
              (do_global_mode ? " (global)" : ""),
              (do_meta_mode ? " (meta)" : ""));
    
    Restrict (cctkGH);
    
    if (reflevel < reflevels-1) {
      ScheduleTraverse (where, "CCTK_POSTRESTRICT", cctkGH);
    }
    
    ScheduleTraverse (where, "CCTK_POSTSTEP", cctkGH);
    
    PoisonCheck (cctkGH, alltimes);
    
    timer.stop();
  }
  
  void
  initialise_3tl_advance_time_2 (cGH * const cctkGH)
  {
    Waypoint ("Advancing time");
    
    cctkGH->cctk_time
      = global_time + 2 * delta_time * mglevelfact / timereflevelfact;
    for (int m=0; m<maps; ++m) {
      vtt.at(m)->advance_time (reflevel, mglevel);
    }
    
    CycleTimeLevels (cctkGH);
  }
  
  void
  initialise_3tl_evolve_Ic (cGH * const cctkGH)
  {
    char const * const where = "Initialise3TL::EvolveIc";
    static Timer timer (where);
    timer.start();
    
    Waypoint ("Initialisation 3TL evolution I (c) (backwards) at iteration"
              " %d time %g%s%s",
              cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
              (do_global_mode ? " (global)" : ""),
              (do_meta_mode ? " (meta)" : ""));
    
    CalculateChecksums (cctkGH, allbutcurrenttime);
    Poison (cctkGH, currenttimebutnotifonly);
    
    // Evolve backward
    ScheduleTraverse (where, "CCTK_PRESTEP", cctkGH);
    ScheduleTraverse (where, "CCTK_EVOL", cctkGH);
    ScheduleTraverse (where, "CCTK_POSTSTEP", cctkGH);
    
    PoisonCheck (cctkGH, alltimes);
    
    timer.stop();
  }
  
  void
  initialise_3tl_reset_time (cGH * const cctkGH)
  {
    Waypoint ("Resetting time");
    
    cctkGH->cctk_time = global_time;
    for (int m=0; m<maps; ++m) {
      vtt.at(m)->set_time (reflevel, mglevel, 0);
    }
  }
  
  
  
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  
  
  
  void
  print_internal_data ()
  {
    DECLARE_CCTK_PARAMETERS;
    
    if (output_internal_data) {
      CCTK_INFO ("Internal data dump:");
      int const oldprecision = cout.precision();
      cout.precision (17);
      cout << "   global_time: " << global_time << endl
           << "   leveltimes: " << leveltimes << endl
           << "   delta_time: " << delta_time << endl;
      cout.precision (oldprecision);
    }
  }
  
  
  
  void ScheduleTraverse (char const * const where, char const * const name,
                         cGH * const cctkGH)
  {
    ostringstream timernamebuf;
    timernamebuf << where << "::" << name;
    string const timername = timernamebuf.str();
    static std::map <string, Timer *> timers;
    Timer * & mapped = timers[timername];
    if (not mapped) {
      mapped = new Timer (timername.c_str());
    }
    Timer & timer = * mapped;
    
    timer.start();
    ostringstream infobuf;
    infobuf << "Scheduling " << name;
    string const info = infobuf.str();
    Checkpoint (info.c_str());
    CCTK_ScheduleTraverse (name, cctkGH, CallFunction);
    timer.stop();
  }
  
  void OutputGH (char const * const where, cGH * const cctkGH)
  {
    ostringstream buf;
    buf << where << "::OutputGH";
    string const timername = buf.str();
    static Timer timer (timername.c_str());
    
    timer.start();
    CCTK_OutputGH (cctkGH);
    timer.stop();
  }
  
  
  
} // namespace Carpet
