#include <cassert>
#include <cstdio>
#include <cstdlib>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Termination.h"

#include "util_String.h"

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
#  include <fcntl.h>
#  include <unistd.h>
#endif

#include "dist.hh"
#include "th.hh"

#include "carpet.hh"
#include "Timers.hh"



namespace Carpet {

  using namespace std;



  static bool do_terminate (const cGH *cgh,
                            const CCTK_REAL time, const int iteration)
  {
    DECLARE_CCTK_PARAMETERS;

    bool term;

    // Early return for non-active reflevels to save the MPI_Allreduce() below
    if (iteration % (maxtimereflevelfact / timereffacts.at(reflevels-1)) != 0)
    {

      return false;

    } else if (terminate_next or CCTK_TerminationReached(cgh)) {

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
                                 and runtime >= 60.0 * max_runtime);
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
      }

    }

    {
      int local, global;
      local = term;
      MPI_Allreduce (&local, &global, 1, MPI_INT, MPI_LOR, dist::comm());
      term = global;
    }

    return term;
  }


  static void AdvanceTime( cGH* cgh, CCTK_REAL initial_time );
  static bool Regrid( cGH* cgh );
  static void PreRegrid( cGH* cgh );
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

    // Set up timers
    static Timer evolve_timer       (timerSet, "Evolve");
    static Timer do_terminate_timer (timerSet, "Evolve::do_terminate");
    static Timer advance_time_timer (timerSet, "Evolve::AdvanceTime");
    static Timer preregrid_timer    (timerSet, "Evolve::PreRegrid");
    static Timer regrid_timer       (timerSet, "Evolve::Regrid");
    static Timer postregrid_timer   (timerSet, "Evolve::PostRegrid");
    static Timer evolution_i_timer  (timerSet, "Evolve::EvolutionI");
    static Timer restrict_timer     (timerSet, "Evolve::Restrict");
    static Timer evolution_ii_timer (timerSet, "Evolve::EvolutionII");

    // Timing statistics
    InitTiming (cgh);

    // Main loop
    evolve_timer.start();
    for (;;) {
      do_terminate_timer.start();
      bool const do_term =
        do_terminate (cgh, cgh->cctk_time, cgh->cctk_iteration);
      do_terminate_timer.stop();
      if (do_term) break;

      advance_time_timer.start();
      AdvanceTime( cgh, cctk_initial_time );
      advance_time_timer.stop();

      if ((cgh->cctk_iteration-1)
          % (maxtimereflevelfact / timereffacts.at(reflevels-1)) == 0) {
        Waypoint ("Evolving iteration %d at t=%g",
                  cgh->cctk_iteration, (double)cgh->cctk_time);
      }

      preregrid_timer.start();
      PreRegrid (cgh);
      preregrid_timer.stop();

      regrid_timer.start();
      bool const did_regrid = Regrid (cgh);
      regrid_timer.stop();
      
      if (did_regrid) {
        postregrid_timer.start();
        PostRegrid (cgh);
        postregrid_timer.stop();
      }

      evolution_i_timer.start();
      EvolutionI (cgh);
      evolution_i_timer.stop();

      restrict_timer.start();
      Evolution_Restrict (cgh);
      restrict_timer.stop();

      evolution_ii_timer.start();
      EvolutionII (cgh);
      evolution_ii_timer.stop();

      if (output_internal_data) {
        CCTK_INFO ("Internal data dump:");
        const int oldprecision = cout.precision();
        cout.precision (17);
        cout << "   global_time: " << global_time << endl
             << "   leveltimes: " << leveltimes << endl
             << "   delta_time: " << delta_time << endl;
        cout.precision (oldprecision);
      }

      // Print timer values
      if ((cgh->cctk_iteration - 1) %
          (maxtimereflevelfact / timereffacts.at(reflevels - 1)) == 0 and
          output_timers_every > 0 and
          (cgh->cctk_iteration - 1) % output_timers_every == 0)
      {
        timerSet.printData (cgh, timer_file);
      }

    } // end main loop
    evolve_timer.stop();

    Waypoint ("Done with evolution loop");

    return 0;
  }

  void AdvanceTime( cGH* cgh, CCTK_REAL initial_time )
  {
    DECLARE_CCTK_PARAMETERS;
    
    ++cgh->cctk_iteration;
    if (! adaptive_stepsize) {
      global_time = initial_time
        + cgh->cctk_iteration * delta_time / maxtimereflevelfact;
      cgh->cctk_time = global_time;
    } else {
      cgh->cctk_time += delta_time;
      global_time = cgh->cctk_time;
    }
  }

  void PreRegrid( cGH* cgh )
  {
    for (int rl=0; rl<reflevels; ++rl) {
      for (int ml=mglevels-1; ml>=0; --ml) {
        // Regridding may change coarser grids, so that postregrid has
        // to be run on all levels.  For symmetry, we also run
        // preregrid on all levels.
        // const int do_every = maxtimereflevelfact / timereffacts.at(rl);
        const int do_every = maxtimereflevelfact / timereffacts.at(reflevels-1);
        if ((cgh->cctk_iteration-1) % do_every == 0) {
          enter_global_mode (cgh, ml);
          enter_level_mode (cgh, rl);

          do_global_mode = reflevel==0;
          do_meta_mode = do_global_mode and mglevel==mglevels-1;

          Waypoint ("Preregrid at iteration %d time %g%s%s",
                    cgh->cctk_iteration, (double)cgh->cctk_time,
                    (do_global_mode ? " (global)" : ""),
                    (do_meta_mode ? " (meta)" : ""));

          Checkpoint ("Scheduling PREREGRID");
          CCTK_ScheduleTraverse ("CCTK_PREREGRID", cgh, CallFunction);

          leave_level_mode (cgh);
          leave_global_mode (cgh);
        }
      }
    }
  }

  bool Regrid( cGH* cgh )
  {
    bool did_regrid = false;

    for (int rl=0; rl<reflevels; ++rl) {
      const int ml=0;
      const int do_every = maxtimereflevelfact / timereffacts.at(rl);
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
        // Regridding may change coarser grids, so that postregrid has
        // to be run on all levels.
        // const int do_every = maxtimereflevelfact / timereffacts.at(rl);
        const int do_every = maxtimereflevelfact / timereffacts.at(reflevels-1);
        if ((cgh->cctk_iteration-1) % do_every == 0) {
          enter_global_mode (cgh, ml);
          enter_level_mode (cgh, rl);

          do_global_mode = reflevel==0;
          do_meta_mode = do_global_mode and mglevel==mglevels-1;

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
  }

  void EvolutionI( cGH* cgh )
  {
    for (int ml=mglevels-1; ml>=0; --ml) {

      bool have_done_global_mode = false;
      bool have_done_anything = false;

      for (int rl=0; rl<reflevels; ++rl) {
        const int do_every
          = ipow(mgfact, ml) * (maxtimereflevelfact / timereffacts.at(rl));
        if ((cgh->cctk_iteration-1) % do_every == 0) {
          enter_global_mode (cgh, ml);
          enter_level_mode (cgh, rl);

          do_global_mode = ! have_done_global_mode;
          do_meta_mode = do_global_mode and mglevel==mglevels-1;
          assert (! (have_done_global_mode and do_global_mode));
          have_done_global_mode |= do_global_mode;
          have_done_anything = true;

          // Advance times
          for (int m=0; m<maps; ++m) {
            vtt.at(m)->advance_time (reflevel, mglevel);
          }
          cgh->cctk_time = (global_time
                            - delta_time / maxtimereflevelfact
                            + delta_time * mglevelfact / timereflevelfact);
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

          // Timing statistics
          StepTiming (cgh);

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
          = ipow(mgfact, ml) * (maxtimereflevelfact / timereffacts.at(rl));
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
    DECLARE_CCTK_PARAMETERS;

    for (int ml=mglevels-1; ml>=0; --ml) {

      bool have_done_global_mode = false;
      bool have_done_anything = false;

      for (int rl=0; rl<reflevels; ++rl) {
        const int do_every
          = ipow(mgfact, ml) * (maxtimereflevelfact / timereffacts.at(rl));
        if (cgh->cctk_iteration % do_every == 0) {
          enter_global_mode (cgh, ml);
          enter_level_mode (cgh, rl);

          int finest_active_reflevel = -1;
          {
            for (int rl_=0; rl_<reflevels; ++rl_) {
              const int do_every_
                = (ipow(mgfact, ml)
                   * (maxtimereflevelfact / timereffacts.at(rl_)));
              if (cgh->cctk_iteration % do_every_ == 0) {
                finest_active_reflevel = rl_;
              }
            }
            assert (finest_active_reflevel >= 0);
          }
          do_global_mode = rl == finest_active_reflevel;
          do_meta_mode = do_global_mode and mglevel==mglevels-1;
          assert (! (have_done_global_mode and do_global_mode));
          have_done_global_mode |= do_global_mode;
          have_done_anything = true;

          Waypoint ("Evolution II at iteration %d time %g%s%s",
                    cgh->cctk_iteration, (double)cgh->cctk_time,
                    (do_global_mode ? " (global)" : ""),
                    (do_meta_mode ? " (meta)" : ""));

          if (rl < reflevels-1) {
            Checkpoint ("Scheduling POSTRESTRICT");
            CCTK_ScheduleTraverse ("CCTK_POSTRESTRICT", cgh, CallFunction);
          }

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

          if (do_global_mode and
              print_timestats_every > 0 and
              cgh->cctk_iteration % print_timestats_every == 0)
          {
            // Timing statistics
            PrintTimingStats (cgh);
          }

          leave_level_mode (cgh);
          leave_global_mode (cgh);
        }
      }

      if (have_done_anything)
        assert (have_done_global_mode);

    }

  }

} // namespace Carpet
