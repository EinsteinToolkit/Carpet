#include <cassert>
#include <cstdlib>
#include <iostream>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctki_GHExtensions.h"
#include "cctki_ScheduleBindings.h"
#include "cctki_WarnLevel.h"

#include "carpet.hh"



namespace Carpet {

  using namespace std;


  static void output_the_grid_structure ( cGH* cgh );
  static void register_coordinates_and_check_parameters ( cGH* cgh );
  static void recovery_I ( cGH* cgh, int rl );
  static void recovery_Regrid ( cGH* cgh, int rl );
  static void recovery_II ( cGH* cgh );
  static void
  initialisation_I ( cGH* cgh, int rl, int init_each_timelevel );
  static void initialise_rewind ( cGH* cgh, int num_tl );
  static void initialise_Schedule_INITIAL ( cGH* cgh, int num_tl );
  static void
  initialise_Regrid ( cGH* cgh, int rl, int prolongate_initial_data );
  static void initialise_Restrict ( cGH* cgh );
  static void initialisation_II ( cGH* cgh );
  static void get_two_extra_timelevels_of_data ( cGH* cgh );
  static void initialise_3_Timelevels ( cGH* cgh );
  static void initialise_Flip_Timelevels ( cGH* cgh );
  static void initialise_evolve_3TL_backwards_Ib ( cGH* cgh );
  static void initialise_evolve_3TL_backwards_IIb_Ic ( cGH* cgh );
  static void initialise_Flip_Timelevels_back ( cGH* cgh );
  static void initialisation_III ( cGH* cgh );

  int Initialise (tFleshConfig* fc)
  {
    DECLARE_CCTK_PARAMETERS;

    const int convlev = 0;
    cGH* const cgh = CCTK_SetupGH (fc, convlev);
    CCTKi_AddGH (fc, convlev, cgh);

    do_global_mode = true;
    do_meta_mode = true;
    global_time = cctk_initial_time;
    delta_time = 1.0;

    cgh->cctk_iteration = 0;
    cgh->cctk_time = global_time;
    cgh->cctk_delta_time = delta_time;

    // Delay checkpoint until MPI has been initialised
    Waypoint ("Starting initialisation");

    CCTKi_ScheduleGHInit (cgh); // Enable storage and communication

    CCTKi_InitGHExtensions (cgh);

    output_the_grid_structure (cgh);

    register_coordinates_and_check_parameters (cgh);

    if (fc->recovered) {

      for (int rl=0; rl<reflevels; ++rl) {
        recovery_I (cgh, rl);
        recovery_Regrid (cgh, rl);
      }

      recovery_II (cgh);

      if (output_internal_data) {
        CCTK_INFO ("Internal data dump:");
        cout << "   global_time: " << global_time << endl
             << "   leveltimes: " << leveltimes << endl
             << "   delta_time: " << delta_time << endl;
      }

    } else {

      for (int rl=0; rl<reflevels; ++rl) {
        initialisation_I (cgh, rl, init_each_timelevel);
        initialise_Regrid (cgh, rl, prolongate_initial_data);
      }

      initialise_Restrict (cgh);

      initialisation_II (cgh);

      if (output_internal_data) {
        CCTK_INFO ("Internal data dump:");
        cout << "   global_time: " << global_time << endl
             << "   leveltimes: " << leveltimes << endl
             << "   delta_time: " << delta_time << endl;
      }

      if (init_3_timelevels) {
        get_two_extra_timelevels_of_data (cgh);
      }
    }

    initialisation_III (cgh);

    if (output_internal_data) {
      CCTK_INFO ("Internal data dump:");
      cout << "   global_time: " << global_time << endl
           << "   leveltimes: " << leveltimes << endl
           << "   delta_time: " << delta_time << endl;
    }

    Waypoint ("Done with initialisation");

    return 0;
  }

  void output_the_grid_structure (cGH* cgh)
  {
    // Loop over maps
    for (int m=0; m<maps; ++m) {
      // Write grid structure to file
      OutputGridStructure
        (cgh, m,
         vhh.at(m)->extents(), vhh.at(m)->outer_boundaries(),
         vhh.at(m)->processors());
    } // loop over maps
  }

  void register_coordinates_and_check_parameters (cGH* cgh)
  {
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
  }

  void recovery_I (cGH* cgh, int rl)
  {
    BEGIN_MGLEVEL_LOOP(cgh) {
      enter_level_mode (cgh, rl);
      do_global_mode = reflevel==0;
      do_meta_mode = do_global_mode and mglevel==mglevels-1;

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
  }

  void recovery_Regrid (cGH* cgh, int rl)
  {
    bool did_regrid = false;
    {
      const int ml=0;
      enter_global_mode (cgh, ml);
      enter_level_mode (cgh, rl);

      // Regrid
      Checkpoint ("Regrid");
      did_regrid |= Regrid (cgh, true, false);

      leave_level_mode (cgh);
      leave_global_mode (cgh);
    } // ml

    if (did_regrid) {
      BEGIN_MGLEVEL_LOOP(cgh) {
        enter_level_mode (cgh, rl);
        do_global_mode = true;
        do_meta_mode = do_global_mode and mglevel==mglevels-1;

        Waypoint ("Postregrid at iteration %d time %g%s%s",
                  cgh->cctk_iteration, (double)cgh->cctk_time,
                  (do_global_mode ? " (global)" : ""),
                  (do_meta_mode ? " (meta)" : ""));

        // Postregrid
        Checkpoint ("Scheduling POSTREGRID");
        CCTK_ScheduleTraverse ("CCTK_POSTREGRID", cgh, CallFunction);

        leave_level_mode (cgh);
      } END_MGLEVEL_LOOP;
    } // if did_regrid
  }

  void recovery_II (cGH* cgh)
  {
    for (int rl=0; rl<reflevels; ++rl) {
      BEGIN_MGLEVEL_LOOP(cgh) {
        enter_level_mode (cgh, rl);
        do_global_mode = reflevel==reflevels-1;
        do_meta_mode = do_global_mode and mglevel==mglevels-1;

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
    }
  }

  void initialisation_I (cGH* cgh, int rl, int init_each_timelevel)
  {
    BEGIN_MGLEVEL_LOOP(cgh) {
      enter_level_mode (cgh, rl);
      do_global_mode = reflevel==0;
      do_meta_mode = do_global_mode and mglevel==mglevels-1;

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

      initialise_rewind (cgh, num_tl);

      initialise_Schedule_INITIAL (cgh, num_tl);

      // Checking
      PoisonCheck (cgh, currenttime);

      leave_level_mode (cgh);
    } END_MGLEVEL_LOOP;
  }

  void initialise_rewind (cGH* cgh, int num_tl)
  {
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
  }

  void initialise_Schedule_INITIAL (cGH* cgh, int num_tl)
  {
    const bool outer_do_global_mode = do_global_mode;
    for (int tl=num_tl-1; tl>=0; --tl) {
      do_global_mode = outer_do_global_mode and tl==0;

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
  }

  void initialise_Regrid (cGH* cgh, int rl, int prolongate_initial_data)
  {
    bool did_regrid = false;
    {
      const int ml=0;
      enter_global_mode (cgh, ml);
      enter_level_mode (cgh, rl);

      // Regrid
      Checkpoint ("Regrid");
      did_regrid |= Regrid (cgh, false, prolongate_initial_data);

      leave_level_mode (cgh);
      leave_global_mode (cgh);
    } // ml

    if (did_regrid) {
      BEGIN_MGLEVEL_LOOP(cgh) {
        enter_level_mode (cgh, rl);
        do_global_mode = true;
        do_meta_mode = do_global_mode and mglevel==mglevels-1;

        Waypoint ("Postregrid at iteration %d time %g%s%s",
                  cgh->cctk_iteration, (double)cgh->cctk_time,
                  (do_global_mode ? " (global)" : ""),
                  (do_meta_mode ? " (meta)" : ""));

        // Postregrid
        Checkpoint ("Scheduling POSTREGRID");
        CCTK_ScheduleTraverse ("CCTK_POSTREGRID", cgh, CallFunction);

        leave_level_mode (cgh);
      } END_MGLEVEL_LOOP;
    } // if did_regrid
  }

  void initialise_Restrict (cGH* cgh)
  {
    for (int rl=reflevels-1; rl>=0; --rl) {
      BEGIN_MGLEVEL_LOOP(cgh) {
        enter_level_mode (cgh, rl);

        Waypoint ("Initialisation/Restrict at iteration %d time %g",
                  cgh->cctk_iteration, (double)cgh->cctk_time);

        Restrict (cgh);

        leave_level_mode (cgh);
      } END_MGLEVEL_LOOP;
    }
  }

  void initialisation_II (cGH* cgh)
  {
    for (int rl=0; rl<reflevels; ++rl) {
      BEGIN_MGLEVEL_LOOP(cgh) {
        enter_level_mode (cgh, rl);
        do_global_mode = reflevel==reflevels-1;
        do_meta_mode = do_global_mode and mglevel==mglevels-1;

        Waypoint ("Initialisation II at iteration %d time %g%s%s",
                  cgh->cctk_iteration, (double)cgh->cctk_time,
                  (do_global_mode ? " (global)" : ""),
                  (do_meta_mode ? " (meta)" : ""));

        if (rl < reflevels-1) {
          Checkpoint ("Scheduling POSTRESTRICTINITIAL");
          CCTK_ScheduleTraverse
            ("CCTK_POSTRESTRICTINITIAL", cgh, CallFunction);
        }

        Checkpoint ("Scheduling POSTINITIAL");
        CCTK_ScheduleTraverse ("CCTK_POSTINITIAL", cgh, CallFunction);

        Checkpoint ("Scheduling POSTSTEP");
        CCTK_ScheduleTraverse ("CCTK_POSTSTEP", cgh, CallFunction);

        PoisonCheck (cgh, alltimes);
        CheckChecksums (cgh, allbutcurrenttime);

        leave_level_mode (cgh);
      } END_MGLEVEL_LOOP;
    }
  }

  // Use Scott Hawley's algorithm for getting two extra timelevels of
  // data
  void get_two_extra_timelevels_of_data (cGH* cgh)
  {
        Waypoint ("Initialising three timelevels");

        initialise_3_Timelevels (cgh);

        delta_time *= -1;
        initialise_Flip_Timelevels (cgh);

        initialise_evolve_3TL_backwards_Ib (cgh);

        Waypoint ("Hourglass structure in place");

        initialise_evolve_3TL_backwards_IIb_Ic (cgh);

        delta_time *= -1;
        initialise_Flip_Timelevels_back (cgh);

        Waypoint ("Finished initialising three timelevels");
  }

  void initialise_3_Timelevels (cGH* cgh)
  {
    for (int rl=0; rl<reflevels; ++rl) {
      BEGIN_MGLEVEL_LOOP(cgh) {
        enter_level_mode (cgh, rl);
        do_global_mode = reflevel==0;
        do_meta_mode = do_global_mode and mglevel==mglevels-1;

        // Advance times
        for (int m=0; m<maps; ++m) {
          vtt.at(m)->advance_time (reflevel, mglevel);
        }
        cgh->cctk_time
          = global_time + delta_time * mglevelfact / reflevelfact;
        CycleTimeLevels (cgh);

        Waypoint ("Initialisation 3TL evolution I (a) (forwards) at iteration"
                  " %d time %g%s%s",
                  cgh->cctk_iteration, (double)cgh->cctk_time,
                  (do_global_mode ? " (global)" : ""),
                  (do_meta_mode ? " (meta)" : ""));

        CalculateChecksums (cgh, allbutcurrenttime);
        Poison (cgh, currenttimebutnotifonly);

        // Evolve forward
        Checkpoint ("Scheduling PRESTEP");
        CCTK_ScheduleTraverse ("CCTK_PRESTEP", cgh, CallFunction);
        Checkpoint ("Scheduling EVOL");
        CCTK_ScheduleTraverse ("CCTK_EVOL", cgh, CallFunction);

        PoisonCheck (cgh, currenttime);

        leave_level_mode (cgh);
      } END_MGLEVEL_LOOP;
    }
  }

  void initialise_Flip_Timelevels (cGH* cgh)
  {
    for (int rl=0; rl<reflevels; ++rl) {
      BEGIN_MGLEVEL_LOOP(cgh) {
        enter_level_mode (cgh, rl);
        do_global_mode = reflevel==0;
        do_meta_mode = do_global_mode and mglevel==mglevels-1;

        // Flip time levels
        Waypoint ("Flipping timelevels");
        FlipTimeLevels (cgh);

        cgh->cctk_time
          = global_time + delta_time * mglevelfact / reflevelfact;

        leave_level_mode (cgh);
      } END_MGLEVEL_LOOP;
    }
  }

  void initialise_evolve_3TL_backwards_Ib (cGH* cgh)
  {
    for (int rl=0; rl<reflevels; ++rl) {
      BEGIN_MGLEVEL_LOOP(cgh) {
        enter_level_mode (cgh, rl);
        do_global_mode = reflevel==0;
        do_meta_mode = do_global_mode and mglevel==mglevels-1;

        Waypoint ("Initialisation 3TL evolution I (b) (backwards) at iteration"
                  " %d time %g%s%s",
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
    }
  }

  // Evolve each level "backwards" one more timestep
  // Starting with the finest level and proceeding to the coarsest
  void initialise_evolve_3TL_backwards_IIb_Ic (cGH* cgh)
  {
    for (int rl=reflevels-1; rl>=0; --rl) {
      BEGIN_MGLEVEL_LOOP(cgh) {
        enter_level_mode (cgh, rl);
        do_global_mode = reflevel==0;
        do_meta_mode = do_global_mode and mglevel==mglevels-1;

        Waypoint ("Initialisation 3TL evolution II (b) (backwards) at iteration"
                  " %d time %g%s%s",
                  cgh->cctk_iteration, (double)cgh->cctk_time,
                  (do_global_mode ? " (global)" : ""),
                  (do_meta_mode ? " (meta)" : ""));

        Restrict (cgh);

        if (rl < reflevels-1) {
          Checkpoint ("Scheduling POSTRESTRICT");
          CCTK_ScheduleTraverse ("CCTK_POSTRESTRICT", cgh, CallFunction);
        }

        Checkpoint ("Scheduling POSTSTEP");
        CCTK_ScheduleTraverse ("CCTK_POSTSTEP", cgh, CallFunction);

        PoisonCheck (cgh, alltimes);

        // Advance times
        for (int m=0; m<maps; ++m) {
          vtt.at(m)->advance_time (reflevel, mglevel);
        }
        cgh->cctk_time
          = global_time + 2 * delta_time * mglevelfact / reflevelfact;
        CycleTimeLevels (cgh);

        Waypoint ("Initialisation 3TL evolution I (c) (backwards) at iteration"
                  " %d time %g%s%s",
                  cgh->cctk_iteration, (double)cgh->cctk_time,
                  (do_global_mode ? " (global)" : ""),
                  (do_meta_mode ? " (meta)" : ""));

        CalculateChecksums (cgh, allbutcurrenttime);
        Poison (cgh, currenttimebutnotifonly);

        // Evolve backward
        Checkpoint ("Scheduling PRESTEP");
        CCTK_ScheduleTraverse ("CCTK_PRESTEP", cgh, CallFunction);
        Checkpoint ("Scheduling EVOL");
        CCTK_ScheduleTraverse ("CCTK_EVOL", cgh, CallFunction);
        Checkpoint ("Scheduling POSTSTEP");
        CCTK_ScheduleTraverse ("CCTK_POSTSTEP", cgh, CallFunction);

        PoisonCheck (cgh, alltimes);

        leave_level_mode (cgh);
      } END_MGLEVEL_LOOP;
    }
  }

  void initialise_Flip_Timelevels_back (cGH* cgh)
  {
    for (int rl=0; rl<reflevels; ++rl) {
      BEGIN_MGLEVEL_LOOP(cgh) {
        enter_level_mode (cgh, rl);
        do_global_mode = reflevel==0;
        do_meta_mode = do_global_mode and mglevel==mglevels-1;

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
    }
  }

  void initialisation_III (cGH* cgh)
  {
    for (int rl=0; rl<reflevels; ++rl) {
      BEGIN_MGLEVEL_LOOP(cgh) {
        enter_level_mode (cgh, rl);
        do_global_mode = reflevel==reflevels-1;
        do_meta_mode = do_global_mode and mglevel==mglevels-1;

        Waypoint ("Initialisation III at iteration %d time %g%s%s",
                  cgh->cctk_iteration, (double)cgh->cctk_time,
                  (do_global_mode ? " (global)" : ""),
                  (do_meta_mode ? " (meta)" : ""));

        const int do_every =
          ipow(mgfact, mglevel) * (maxreflevelfact / ipow(reffact, rl));
        if (cgh->cctk_iteration % do_every == 0)
        {
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
        }

        leave_level_mode (cgh);
      } END_MGLEVEL_LOOP;
    }
  }

} // namespace Carpet
