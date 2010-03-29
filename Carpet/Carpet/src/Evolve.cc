#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <sstream>

#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctk_Termination.h>

#include <util_String.h>

#include <dist.hh>
#include <th.hh>

#include "carpet.hh"
#include "Timers.hh"



namespace Carpet {
  
  using namespace std;
  
  
  
  static bool do_terminate (const cGH * cctkGH);
  
  static void AdvanceTime (cGH * cctkGH);
  static void CallRegrid (cGH * cctkGH);
  static void CallEvol (cGH * cctkGH);
  static void CallRestrict (cGH * cctkGH);
  static void CallAnalysis (cGH * cctkGH);
  
  static void print_internal_data ();
  
  static void ScheduleTraverse
  (char const * where, char const * name, cGH * cctkGH);
  static void OutputGH (char const * where, cGH * cctkGH);
  
  
  
  int
  Evolve (tFleshConfig * const fc)
  {
    DECLARE_CCTK_PARAMETERS;
    
    Waypoint ("Starting evolution loop");
    
    int const convlev = 0;
    cGH* cctkGH = fc->GH[convlev];
    
    // Tapered grids
    do_taper = use_tapered_grids;
    
    // Main loop
    BeginTiming (cctkGH);
    static Timer timer ("Evolve");
    timer.start();
    while (not do_terminate (cctkGH)) {
      
      AdvanceTime (cctkGH);
      {
        int const do_every = maxtimereflevelfact / timereffacts.at(reflevels-1);
        if ((cctkGH->cctk_iteration - 1) % do_every == 0) {
          bool const old_do_taper = do_taper;
          do_taper = false;
          ENTER_GLOBAL_MODE (cctkGH, 0) {
            BEGIN_REFLEVEL_LOOP (cctkGH) {
              CallRegrid (cctkGH);
            } END_REFLEVEL_LOOP;
          } LEAVE_GLOBAL_MODE;
          do_taper = old_do_taper;
        }
      }
      CallEvol (cctkGH);
      CallRestrict (cctkGH);
      CallAnalysis (cctkGH);
      print_internal_data ();
      
      // Print timer values
      {
        int const do_every = maxtimereflevelfact / timereffacts.at(reflevels-1);
        if (output_timers_every > 0 and
            cctkGH->cctk_iteration % output_timers_every == 0 and
            cctkGH->cctk_iteration % do_every == 0)
        {
          TimerSet::writeData (cctkGH, timer_file);
        }
      }
      
      // Ensure that all levels have consistent times
      {
        CCTK_REAL const eps = 1.0e-12;
        assert (abs (cctkGH->cctk_time - global_time) <= eps * global_time);
        for (int ml=0; ml<mglevels; ++ml) {
          for (int rl=0; rl<reflevels; ++rl) {
            int const do_every =
              ipow (mgfact, ml) * (maxtimereflevelfact / timereffacts.at(rl));
            if (cctkGH->cctk_iteration % do_every == 0) {
              assert (abs (leveltimes.at(ml).at(rl) - global_time) <=
                      eps * global_time);
            }
          }
        }
      }
      
    } // end main loop
    timer.stop();
    
    do_taper = false;
    
    Waypoint ("Done with evolution loop");
    
    return 0;
  }
  
  
  
  bool
  do_terminate (const cGH *cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    static Timer timer ("Evolve::do_terminate");
    timer.start();
    
    bool term;
    
    // Do not test on non-active reflevels to save the call to
    // MPI_Allreduce below
    int const do_every = maxtimereflevelfact / timereffacts.at(reflevels-1);
    if (cctkGH->cctk_iteration % do_every != 0)
    {
      
      term = false;
      
    } else {
      
      if (terminate_next or CCTK_TerminationReached(cctkGH)) {
        
        // Terminate if someone or something said so
        term = true;
        
      } else {
        
        // Test the various conditions
        bool const term_iter
          = cctk_itlast >= 0 and cctkGH->cctk_iteration >= cctk_itlast;
        bool const term_time
          = (delta_time > 0.0
             ? (cctkGH->cctk_time
                >= cctk_final_time - 1.0e-8 * cctkGH->cctk_delta_time)
             : (cctkGH->cctk_time
                <= cctk_final_time + 1.0e-8 * cctkGH->cctk_delta_time));
        bool const term_runtime
          = max_runtime > 0.0 and CCTK_RunTime() >= 60.0 * max_runtime;
        
        if (CCTK_Equals(terminate, "never")) {
          term = false;
        } else if (CCTK_Equals(terminate, "iteration")) {
          term = term_iter;
        } else if (CCTK_Equals(terminate, "time")) {
          term = term_time;
        } else if (CCTK_Equals(terminate, "runtime")) {
          term = term_runtime;
        } else if (CCTK_Equals(terminate, "any")) {
          term = term_iter or term_time or term_runtime;
        } else if (CCTK_Equals(terminate, "all")) {
          term = term_iter and term_time and term_runtime;
        } else if (CCTK_Equals(terminate, "either")) {
          term = term_iter or term_time;
        } else if (CCTK_Equals(terminate, "both")) {
          term = term_iter and term_time;
        } else if (CCTK_Equals(terminate, "immediately")) {
          term = true;
        } else {
          CCTK_WARN (0, "Unsupported termination condition");
          abort ();             // keep the compiler happy
        }
        
      }
      
      // Reduce termination condition
      int local, global;
      local = term;
      MPI_Allreduce (&local, &global, 1, MPI_INT, MPI_LOR, dist::comm());
      term = global;
      
    }
    
    timer.stop();
    return term;
  }
  
  
  
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  
  
  
  void
  AdvanceTime (cGH * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    static Timer timer ("Evolve::AdvanceTime");
    timer.start();
    
    Checkpoint ("AdvanceTime");
    
    ++ cctkGH->cctk_iteration;
    
    if (not adaptive_stepsize) {
      // Avoid accumulation of errors
      global_time = cctk_initial_time
        + cctkGH->cctk_iteration * delta_time / maxtimereflevelfact;
      cctkGH->cctk_time = global_time;
    } else {
      // Take varying step sizes into account
      cctkGH->cctk_time += delta_time;
      global_time = cctkGH->cctk_time;
    }
    
    if ((cctkGH->cctk_iteration-1)
        % (maxtimereflevelfact / timereffacts.at(reflevels-1)) == 0) {
      Waypoint ("Evolving iteration %d at t=%g",
                cctkGH->cctk_iteration, (double)cctkGH->cctk_time);
    }
    
    timer.stop();
  }
  
  
  
  void
  CallRegrid (cGH * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    char const * const where = "Evolve::CallRegrid";
    static Timer timer (where);
    timer.start();
    
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
    
    Waypoint ("Preregrid at iteration %d time %g%s%s",
              cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
              (do_global_mode ? " (global)" : ""),
              (do_meta_mode ? " (meta)" : ""));
    
    // Preregrid
    ScheduleTraverse (where, "CCTK_PREREGRID", cctkGH);
    
    // Regrid
    Checkpoint ("Regrid");
    int const oldreflevels = reflevels;
    bool const did_regrid = Regrid (cctkGH, false);
    bool const did_remove_level = reflevels < oldreflevels;
    assert (not did_remove_level or did_regrid);
    
    if (did_regrid) {
      BEGIN_META_MODE (cctkGH) {
        for (int rl=0; rl<reflevels; ++rl) {
          
          bool const did_recompose = Recompose (cctkGH, rl, true);
          
          // Do not omit the global mode call when the finest level
          // does not change:
          // if (did_recompose or (did_remove_level and rl == reflevels - 1)) {
          if (did_recompose or rl == reflevels - 1) {
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
                
                int const num_tl = prolongation_order_time+1;
                
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
  CallEvol (cGH * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    char const * const where = "Evolve::CallEvol";
    static Timer timer (where);
    timer.start();
    
    for (int ml=mglevels-1; ml>=0; --ml) {
      
      bool have_done_global_mode = false;
      bool have_done_anything = false;
      
      for (int rl=0; rl<reflevels; ++rl) {
        int const do_every
          = ipow(mgfact, ml) * (maxtimereflevelfact / timereffacts.at(rl));
        if ((cctkGH->cctk_iteration-1) % do_every == 0) {
          ENTER_GLOBAL_MODE (cctkGH, ml) {
            ENTER_LEVEL_MODE (cctkGH, rl) {
              
              do_early_global_mode = not have_done_global_mode;
              do_late_global_mode = reflevel==reflevels-1;
              do_early_meta_mode = do_early_global_mode and mglevel==mglevels-1;
              do_late_meta_mode = do_late_global_mode and mglevel==0;
              do_global_mode = do_early_global_mode;
              do_meta_mode = do_early_meta_mode;
              assert (not (have_done_global_mode and do_global_mode));
              have_done_global_mode |= do_global_mode;
              have_done_anything = true;
              
              // Advance times
              cctkGH->cctk_time
                = (global_time
                   - delta_time / maxtimereflevelfact
                   + delta_time * mglevelfact / timereflevelfact);
              CCTK_REAL const carpet_time = cctkGH->cctk_time / delta_time;
              for (int m=0; m<maps; ++m) {
                vtt.at(m)->advance_time (reflevel, mglevel);
                if (not adaptive_stepsize) {
#if 0
                  // We must not perform this check, since the
                  // relative accuracy of incrementally adding to the
                  // current time cannot be good enough.  Just setting
                  // the time (see below) is fine.
                  CCTK_REAL const eps = 1.0e-12;
                  static_assert (abs(0.1) > 0,
                                 "Function CarpetLib::abs has wrong signature");
                  CCTK_REAL const level_time =
                    vtt.at(m)->get_time (reflevel, mglevel);
                  if (not (abs (level_time - carpet_time) <=
                           eps * max (carpet_time, 1.0))) {
                    int const oldprecision = cerr.precision();
                    cerr.precision (17);
                    cerr << "ml: " << ml << endl
                         << "rl: " << rl << endl
                         << "m: " << m << endl
                         << "level_time: " << level_time << endl
                         << "carpet_time: " << carpet_time << endl
                         << "(level_time - carpet_time): " << (level_time - carpet_time) << endl;
                    cerr.precision (oldprecision);
                  }
                  assert (abs (level_time - carpet_time) <=
                          eps * max (carpet_time, 1.0));
#endif
                  vtt.at(m)->set_time (reflevel, mglevel, carpet_time);
                }
              }
              CycleTimeLevels (cctkGH);
              
              Waypoint ("Evolution I at iteration %d time %g%s%s",
                        cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
                        (do_global_mode ? " (global)" : ""),
                        (do_meta_mode ? " (meta)" : ""));
              
              // Checking
              CalculateChecksums (cctkGH, allbutcurrenttime);
              Poison (cctkGH, currenttimebutnotifonly);
              
              // Evolve
              ScheduleTraverse (where, "CCTK_PRESTEP", cctkGH);
              ScheduleTraverse (where, "CCTK_EVOL", cctkGH);
              
              // Checking
              PoisonCheck (cctkGH, currenttime);
              
              // Timing statistics
              StepTiming (cctkGH);
              
            } LEAVE_LEVEL_MODE;
          } LEAVE_GLOBAL_MODE;
        } // if do_every
      }   // for rl
      
      if (have_done_anything) assert (have_done_global_mode);
      
    } // for ml
    
    timer.stop();
  }
  
  
  
  void
  CallRestrict (cGH * const cctkGH)
  {
    static Timer timer ("Evolve::CallRestrict");
    timer.start();
    
    for (int ml=mglevels-1; ml>=0; --ml) {
      for (int rl=reflevels-1; rl>=0; --rl) {
        int const do_every
          = ipow(mgfact, ml) * (maxtimereflevelfact / timereffacts.at(rl));
        if (cctkGH->cctk_iteration % do_every == 0) {
          ENTER_GLOBAL_MODE (cctkGH, ml) {
            ENTER_LEVEL_MODE (cctkGH, rl) {
              
              Waypoint ("Evolution/Restrict at iteration %d time %g",
                        cctkGH->cctk_iteration, (double)cctkGH->cctk_time);
              
              Restrict (cctkGH);
              
            } LEAVE_LEVEL_MODE;
          } LEAVE_GLOBAL_MODE;
        } // if do_every
      }   // for rl
    }     // for ml
    
    timer.stop();
  }
  
  
  
  void
  CallAnalysis  (cGH * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    char const * const where = "Evolve::CallAnalysis";
    static Timer timer (where);
    timer.start();
    
    for (int ml=mglevels-1; ml>=0; --ml) {
      
      bool have_done_global_mode = false;
      bool have_done_early_global_mode = false;
      bool have_done_late_global_mode = false;
      bool have_done_anything = false;
      
      for (int rl=0; rl<reflevels; ++rl) {
        int const do_every
          = ipow(mgfact, ml) * (maxtimereflevelfact / timereffacts.at(rl));
        if (cctkGH->cctk_iteration % do_every == 0) {
          ENTER_GLOBAL_MODE (cctkGH, ml) {
            ENTER_LEVEL_MODE (cctkGH, rl) {
              
              do_early_global_mode = not have_done_early_global_mode;
              do_late_global_mode = reflevel==reflevels-1;
              do_early_meta_mode = do_early_global_mode and mglevel==mglevels-1;
              do_late_meta_mode = do_late_global_mode and mglevel==0;
              do_global_mode = do_late_global_mode;
              do_meta_mode = do_global_mode and do_late_meta_mode;
              assert (not (have_done_global_mode and do_global_mode));
              assert (not (have_done_early_global_mode and
                           do_early_global_mode));
              assert (not (have_done_late_global_mode and
                           do_late_global_mode));
              have_done_global_mode |= do_global_mode;
              have_done_early_global_mode |= do_early_global_mode;
              have_done_late_global_mode |= do_late_global_mode;
              have_done_anything = true;
              
              Waypoint ("Evolution II at iteration %d time %g%s%s",
                        cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
                        (do_global_mode ? " (global)" : ""),
                        (do_meta_mode ? " (meta)" : ""));
              
              if (reflevel < reflevels-1) {
                ScheduleTraverse (where, "CCTK_POSTRESTRICT", cctkGH);
              }
              
              // Poststep
              ScheduleTraverse (where, "CCTK_POSTSTEP", cctkGH);
              
              // Checking
              PoisonCheck (cctkGH, currenttime);
              CalculateChecksums (cctkGH, currenttime);
              
              // Checkpoint
              ScheduleTraverse (where, "CCTK_CHECKPOINT", cctkGH);
              
              // Analysis
              ScheduleTraverse (where, "CCTK_ANALYSIS", cctkGH);
              
              if (do_late_global_mode) {
                // Timing statistics
                UpdateTimingStats (cctkGH);
              }
              
              // Output
              OutputGH (where, cctkGH);
              
              // Checking
              CheckChecksums (cctkGH, alltimes);
              
              if (do_late_global_mode) {
                // Timing statistics
                PrintTimingStats (cctkGH);
              }
              
            } LEAVE_LEVEL_MODE;
          } LEAVE_GLOBAL_MODE;
        } // if do_every
      }   // for rl
      
      if (have_done_anything) assert (have_done_global_mode);
      if (have_done_anything) assert (have_done_early_global_mode);
      if (have_done_anything) assert (have_done_late_global_mode);
      
    } // for ml
    
    timer.stop();
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
