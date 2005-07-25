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



  static void output_grid_structure (cGH * cgh);
  static void setup_I (cGH * cgh);

  static void recover_I (cGH * cgh, int rl);
  static void recover_Regrid (cGH * cgh, int rl);
  static void recover_II (cGH * cgh, int rl);

  static void initialise_I (cGH * cgh, int rl);
  static void initialise_I_rewind (cGH * cgh, int num_tl);
  static void initialise_I_initialise (cGH * cgh, int num_tl);
  static void initialise_Regrid (cGH * cgh, int rl);
  static void initialise_Restrict (cGH * cgh, int rl);
  static void initialise_II (cGH * cgh, int rl);
  static void initialise_III (cGH * cgh, int rl);

  static void initialise_3tl (cGH * cgh);
  static void initialise_3tl_advance_time (cGH * cgh);
  static void initialise_3tl_evolve_Ia (cGH * cgh);
  static void initialise_3tl_flip_timelevels (cGH * cgh);
  static void initialise_3tl_evolve_Ib (cGH * cgh);
  static void initialise_3tl_evolve_IIb (cGH * cgh);
  static void initialise_3tl_advance_time_2 (cGH * cgh);
  static void initialise_3tl_evolve_Ic (cGH * cgh);
  static void initialise_3tl_flip_timelevels_back (cGH * cgh);

  static void print_internal_data ();



  int Initialise (tFleshConfig * const fc)
  {
    DECLARE_CCTK_PARAMETERS;

    const int convlev = 0;
    cGH * const cgh = CCTK_SetupGH (fc, convlev);
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

    output_grid_structure (cgh);

    setup_I (cgh);

    if (fc->recovered) {
      // Read data from a checkpoint file

      for (int rl=0; rl<reflevels; ++rl) {
        recover_I (cgh, rl);
        recover_Regrid (cgh, rl);
      }

      for (int rl=0; rl<reflevels; ++rl) {
        recover_II (cgh, rl);
      }

      print_internal_data ();

    } else {
      // Calculate initial data

      for (int rl=0; rl<reflevels; ++rl) {
        initialise_I (cgh, rl);
        initialise_Regrid (cgh, rl);
      }

      for (int rl=reflevels-1; rl>=0; --rl) {
        initialise_Restrict (cgh, rl);
      }

      for (int rl=0; rl<reflevels; ++rl) {
        initialise_II (cgh, rl);
      }

      print_internal_data ();

      if (init_3_timelevels) {
#warning "TODO: ensure that there are 3 timelevels"
#warning "TODO: ensure that prolongation_order_time = 2"
        initialise_3tl (cgh);
      }
    }

    // Analyse initial data
    for (int rl=0; rl<reflevels; ++rl) {
      initialise_III (cgh, rl);
    }

    print_internal_data ();

    Waypoint ("Done with initialisation");

    return 0;
  }



  void output_grid_structure (cGH * const cgh)
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

  void setup_I (cGH * const cgh)
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
    } END_MGLEVEL_LOOP;

    CCTKi_FinaliseParamWarn();
  }



  void recover_I (cGH * const cgh, int const rl)
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

  void recover_Regrid (cGH * const cgh, int const rl)
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

  void recover_II (cGH * const cgh, int const rl)
  {
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



  void initialise_I (cGH * const cgh, int const rl)
  {
    DECLARE_CCTK_PARAMETERS;
    
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

      initialise_I_rewind (cgh, num_tl);

      initialise_I_initialise (cgh, num_tl);

      // Checking
      PoisonCheck (cgh, currenttime);

      leave_level_mode (cgh);
    } END_MGLEVEL_LOOP;
  }

  void initialise_I_rewind (cGH * const cgh, int const num_tl)
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

  void initialise_I_initialise (cGH * const cgh, int const num_tl)
  {
    const bool outer_do_global_mode = do_global_mode;
    for (int tl=num_tl-1; tl>=0; --tl) {
      do_global_mode = outer_do_global_mode and tl==0;

      // Advance times
      for (int m=0; m<maps; ++m) {
        vtt.at(m)->advance_time (reflevel, mglevel);
      }
      cgh->cctk_time
        = global_time - tl * delta_time * mglevelfact / timereflevelfact;
      CycleTimeLevels (cgh);

      // Set up the initial data
      Checkpoint ("Scheduling INITIAL");
      CCTK_ScheduleTraverse ("CCTK_INITIAL", cgh, CallFunction);

    } // for tl
    do_global_mode = outer_do_global_mode;
  }

  void initialise_Regrid (cGH * const cgh, int const rl)
  {
    DECLARE_CCTK_PARAMETERS;
    
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

  void initialise_Restrict (cGH * const cgh, int const rl)
  {
    BEGIN_MGLEVEL_LOOP(cgh) {
      enter_level_mode (cgh, rl);
      
      Waypoint ("Initialisation/Restrict at iteration %d time %g",
                cgh->cctk_iteration, (double)cgh->cctk_time);
      
      Restrict (cgh);
      
      leave_level_mode (cgh);
    } END_MGLEVEL_LOOP;
  }

  void initialise_II (cGH * const cgh, int const rl)
  {
    BEGIN_MGLEVEL_LOOP(cgh) {
      enter_level_mode (cgh, rl);
      do_global_mode = reflevel==reflevels-1;
      do_meta_mode = do_global_mode and mglevel==mglevels-1;
      
      Waypoint ("Initialisation II at iteration %d time %g%s%s",
                cgh->cctk_iteration, (double)cgh->cctk_time,
                (do_global_mode ? " (global)" : ""),
                (do_meta_mode ? " (meta)" : ""));
      
      if (reflevel < reflevels-1) {
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

  void initialise_III (cGH * const cgh, int const rl)
  {
    BEGIN_MGLEVEL_LOOP(cgh) {
      enter_level_mode (cgh, rl);
      do_global_mode = reflevel==reflevels-1;
      do_meta_mode = do_global_mode and mglevel==mglevels-1;
      
      Waypoint ("Initialisation III at iteration %d time %g%s%s",
                cgh->cctk_iteration, (double)cgh->cctk_time,
                (do_global_mode ? " (global)" : ""),
                (do_meta_mode ? " (meta)" : ""));
      
      const int do_every =
        ipow(mgfact, mglevel) * (maxtimereflevelfact / timereffacts.at(rl));
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



  // Use Scott Hawley's algorithm to get two extra timelevels of data
  void initialise_3tl (cGH * const cgh)
  {
    Waypoint ("Initialising three timelevels");

    for (int rl=0; rl<reflevels; ++rl) {
      BEGIN_MGLEVEL_LOOP(cgh) {
        enter_level_mode (cgh, rl);
        do_global_mode = reflevel==0;
        do_meta_mode = do_global_mode and mglevel==mglevels-1;

        initialise_3tl_advance_time (cgh);
        initialise_3tl_evolve_Ia (cgh);

        leave_level_mode (cgh);
      } END_MGLEVEL_LOOP;
    }

    initialise_3tl_flip_timelevels (cgh);

    for (int rl=0; rl<reflevels; ++rl) {
      BEGIN_MGLEVEL_LOOP(cgh) {
        enter_level_mode (cgh, rl);
        do_global_mode = reflevel==0;
        do_meta_mode = do_global_mode and mglevel==mglevels-1;

        initialise_3tl_evolve_Ib (cgh);

        leave_level_mode (cgh);
      } END_MGLEVEL_LOOP;
    }

    Waypoint ("Hourglass structure in place");

    for (int rl=reflevels-1; rl>=0; --rl) {
      BEGIN_MGLEVEL_LOOP(cgh) {
        enter_level_mode (cgh, rl);
        do_global_mode = reflevel==0;
        do_meta_mode = do_global_mode and mglevel==mglevels-1;

        initialise_3tl_evolve_IIb (cgh);
        initialise_3tl_advance_time_2 (cgh);
        initialise_3tl_evolve_Ic (cgh);

        leave_level_mode (cgh);
      } END_MGLEVEL_LOOP;
    }

    initialise_3tl_flip_timelevels_back (cgh);

    Waypoint ("Finished initialising three timelevels");
  }

  void initialise_3tl_advance_time (cGH * const cgh)
  {
    Waypoint ("Advancing time");
    
    for (int m=0; m<maps; ++m) {
      vtt.at(m)->advance_time (reflevel, mglevel);
    }
    cgh->cctk_time = global_time + delta_time * mglevelfact / timereflevelfact;
    
    CycleTimeLevels (cgh);
  }

  void initialise_3tl_evolve_Ia (cGH * const cgh)
  {
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
  }

  void initialise_3tl_flip_timelevels (cGH * const cgh)
  {
    Waypoint ("Flipping timelevels");

    delta_time *= -1;

    BEGIN_MGLEVEL_LOOP(cgh) {
      BEGIN_REFLEVEL_LOOP (cgh) {

        for (int m=0; m<maps; ++m) {
          vtt.at(m)->set_delta
            (reflevel, mglevel, - vtt.at(m)->get_delta (reflevel, mglevel));
          vtt.at(m)->advance_time (reflevel, mglevel);
          vtt.at(m)->advance_time (reflevel, mglevel);
        }
        cgh->cctk_time
          = global_time + delta_time * mglevelfact / timereflevelfact;

        FlipTimeLevels (cgh);

      } END_REFLEVEL_LOOP;
    } END_MGLEVEL_LOOP;
  }

  void initialise_3tl_evolve_Ib (cGH * const cgh)
  {
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
  }

  // Evolve backwards one more timestep
  // Starting with the finest level and proceeding to the coarsest
  void initialise_3tl_evolve_IIb (cGH * const cgh)
  {
    Waypoint ("Initialisation 3TL evolution II (b) (backwards) at iteration"
              " %d time %g%s%s",
              cgh->cctk_iteration, (double)cgh->cctk_time,
              (do_global_mode ? " (global)" : ""),
              (do_meta_mode ? " (meta)" : ""));

    Restrict (cgh);

    if (reflevel < reflevels-1) {
      Checkpoint ("Scheduling POSTRESTRICT");
      CCTK_ScheduleTraverse ("CCTK_POSTRESTRICT", cgh, CallFunction);
    }

    Checkpoint ("Scheduling POSTSTEP");
    CCTK_ScheduleTraverse ("CCTK_POSTSTEP", cgh, CallFunction);

    PoisonCheck (cgh, alltimes);
  }

  void initialise_3tl_advance_time_2 (cGH * const cgh)
  {
    Waypoint ("Advancing time");
    
    for (int m=0; m<maps; ++m) {
      vtt.at(m)->advance_time (reflevel, mglevel);
    }
    cgh->cctk_time
      = global_time + 2 * delta_time * mglevelfact / timereflevelfact;
    
    CycleTimeLevels (cgh);
  }

  void initialise_3tl_evolve_Ic (cGH * const cgh)
  {
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
  }
  
  void initialise_3tl_flip_timelevels_back (cGH * const cgh)
  {
    Waypoint ("Flipping timelevels back");

    delta_time *= -1;

    BEGIN_MGLEVEL_LOOP(cgh) {
      BEGIN_REFLEVEL_LOOP (cgh) {

        for (int m=0; m<maps; ++m) {
          vtt.at(m)->set_delta
            (reflevel, mglevel, - vtt.at(m)->get_delta (reflevel, mglevel));
          vtt.at(m)->advance_time (reflevel, mglevel);
          vtt.at(m)->advance_time (reflevel, mglevel);
        }
        cgh->cctk_time = global_time;

        FlipTimeLevels (cgh);

      } END_REFLEVEL_LOOP;
    } END_MGLEVEL_LOOP;
  }



  void print_internal_data ()
  {
    DECLARE_CCTK_PARAMETERS;

    if (output_internal_data) {
      CCTK_INFO ("Internal data dump:");
      const int oldprecision = cout.precision();
      cout.precision (17);
      cout << "   global_time: " << global_time << endl
           << "   leveltimes: " << leveltimes << endl
           << "   delta_time: " << delta_time << endl;
      cout.precision (oldprecision);
    }
  }

} // namespace Carpet
