#include <assert.h>
#include <stdlib.h>

#include <iostream>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctki_GHExtensions.h"
#include "cctki_ScheduleBindings.h"
#include "cctki_WarnLevel.h"

#include "carpet.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Initialise.cc,v 1.19 2002/09/25 19:55:05 schnetter Exp $";

CCTK_FILEVERSION(Carpet_Initialise_cc)



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
    Waypoint ("starting Initialise...");
    
    // Initialise stuff
    cgh->cctk_iteration = 0;
    cgh->cctk_time = cctk_initial_time;
    
    // Enable storage and communtication
    CCTKi_ScheduleGHInit (cgh);
    
    // Initialise stuff
    CCTKi_InitGHExtensions (cgh);
    
    // Check parameters
    set_mglevel (cgh, 0);
    Waypoint ("Current time is %g", cgh->cctk_time);
    Waypoint ("PARAMCHECK");
    CCTK_ScheduleTraverse ("CCTK_PARAMCHECK", cgh, CallFunction);
    CCTKi_FinaliseParamWarn();
    set_mglevel (cgh, -1);
    
    Waypoint ("Initialising iteration %d...", cgh->cctk_iteration);
    
    BEGIN_REFLEVEL_LOOP(cgh) {
      
      BEGIN_MGLEVEL_LOOP(cgh) {
	
	Waypoint ("%*sCurrent time is %g", 2*reflevel, "", cgh->cctk_time);
	
	// Checking
	Poison (cgh, alltimes);
	
	// Set up the grid
	Waypoint ("%*sScheduling BASEGRID", 2*reflevel, "");
	CCTK_ScheduleTraverse ("CCTK_BASEGRID", cgh, CallFunction);
	
	// Allow the time step to be changed
	if (reflevel==0) {
	  // Initialise time and time step on coarse grid
	  base_delta_time = cgh->cctk_delta_time;
	} else {
// 	assert (abs(cgh->cctk_delta_time - base_delta_time / reflevelfactor)
// 		< 1e-6 * base_delta_time);
	  // This circumvents a bug in CactusBase/Time
	  cgh->cctk_delta_time = base_delta_time / reflevelfact;
	}
	
	// Set up the initial data
	Waypoint ("%*sScheduling INITIAL", 2*reflevel, "");
	CCTK_ScheduleTraverse ("CCTK_INITIAL", cgh, CallFunction);
	
      } END_MGLEVEL_LOOP(cgh);
      
      // Regrid
      Waypoint ("%*sRegrid", 2*reflevel, "");
      Regrid (cgh);
      
    } END_REFLEVEL_LOOP(cgh);
    
    
    
    // Scott's algorithm for getting two extra timelevels of data
    if (init_3_timelevels) {
      int time_dir = 1; //Positive = forward (+t), Negative = backward (-t)
      CCTK_INFO("Initializing THREE timelevels");
      //SHH:check the current time
      BEGIN_REFLEVEL_LOOP(cgh) {
        BEGIN_MGLEVEL_LOOP(cgh) {
	  cout << "Time at rl=" << reflevel << ", ml=" << mglevel << ", is "
	       << tt->get_time (reflevel, mglevel) << endl;
        } END_MGLEVEL_LOOP(cgh);
      } END_REFLEVEL_LOOP(cgh);

      // Stupid hack: Copy data at current time to all prev levels
      BEGIN_REFLEVEL_LOOP(cgh) {
	BEGIN_MGLEVEL_LOOP(cgh) {
	  CopyCurrToPrevTimeLevels(cgh,reflevel);
	} END_MGLEVEL_LOOP(cgh);
      } END_REFLEVEL_LOOP(cgh);
      
      BEGIN_REFLEVEL_LOOP(cgh) {
	BEGIN_MGLEVEL_LOOP(cgh) {
	  // Evolve "forward" (which may be backward for lev=1,3,5,7...)

	  cout << "`forward': Before CycleTimeLevels & AdvanceTime, "
	       << "time on rl=" << reflevel << ", ml=" << mglevel << ", is "
	       << tt->get_time (reflevel, mglevel) << endl;
	  // Cycle time levels
	  // erik: what about arrays?
	  CycleTimeLevels (cgh);
	   
	  // Advance level times
	  tt->advance_time (reflevel, mglevel);
	  cout << "`forward': After CycleTimeLevels & AdvanceTime, "
	       << "time on rl=" << reflevel << ", ml=" << mglevel << ", is "
	       << tt->get_time (reflevel, mglevel) << endl;
	  Waypoint ("%*sScheduling PRESTEP", 2*reflevel, "");
	  CCTK_ScheduleTraverse ("CCTK_PRESTEP", cgh, CallFunction);
	  Waypoint ("%*sScheduling EVOL", 2*reflevel, "");
	  CCTK_ScheduleTraverse ("CCTK_EVOL", cgh, CallFunction);
	  Waypoint ("%*sScheduling POSTSTEP", 2*reflevel, "");
	  CCTK_ScheduleTraverse ("CCTK_POSTSTEP", cgh, CallFunction);
	  
	  // erik: what about arrays?
	  //SHH: There was an error in the way I wrote it.
	  // FlipTimeLevels should ONLY flip on levels reflevel and lower
	  //FlipTimeLevelsOnCoarser(cgh,reflevel);
	  FlipTimeLevels(cgh);
	  for (int rl=0; rl<hh->reflevels(); ++rl) {
	      tt->set_time (rl, mglevel, -1*tt->get_time (rl, mglevel) );
	      tt->set_delta (rl, mglevel,-1*tt->get_delta (rl, mglevel));
	  }
	  // Keep track of which direction (in time) we're integrating
	  time_dir *= -1;
	  cgh->cctk_delta_time = time_dir * base_delta_time / reflevelfact * mglevelfact;
	  
	  cout << "`backward': Before CycleTimeLevels & AdvanceTime, "
	       << "time on rl=" << reflevel << ", ml=" << mglevel << ", is "
	       << tt->get_time (reflevel, mglevel) << endl;
	  // Evolve in the opposite time-direction
	  // Cycle time levels
	  // erik: what about arrays?
//	  CycleTimeLevels (cgh);
	  // Advance level times
//	  tt->advance_time (reflevel, mglevel);
	  cout << "`backward': After CycleTimeLevels & AdvanceTime, "
	       << "time on rl=" << reflevel << ", ml=" << mglevel << ", is "
	       << tt->get_time (reflevel, mglevel) << endl;
	  Waypoint ("%*sScheduling PRESTEP", 2*reflevel, "");
	  CCTK_ScheduleTraverse ("CCTK_PRESTEP", cgh, CallFunction);
	  Waypoint ("%*sScheduling EVOL", 2*reflevel, "");
	  CCTK_ScheduleTraverse ("CCTK_EVOL", cgh, CallFunction);
	  Waypoint ("%*sScheduling POSTSTEP", 2*reflevel, "");
	  CCTK_ScheduleTraverse ("CCTK_POSTSTEP", cgh, CallFunction);
	  
	} END_MGLEVEL_LOOP(cgh);
      } END_REFLEVEL_LOOP(cgh);
      
      //SHH:check the current time
      BEGIN_REFLEVEL_LOOP(cgh) {
        BEGIN_MGLEVEL_LOOP(cgh) {
	  cout << "Time at rl=" << reflevel << ", ml=" << mglevel << ", is "
	       << tt->get_time (reflevel, mglevel) << endl;
        } END_MGLEVEL_LOOP(cgh);
      } END_REFLEVEL_LOOP(cgh);


      
      // Make sure we're pointed backwards, in order to get 2 "previous"
      //   timelevels.  We could change the if statement to 
      //   "if (mod(MaxLevels,2) == 0)", but I prefer to check time_dir 
      //   explicitly, because it's easier to follow and I don't have to
      //   worry about having made a mistake
      if (time_dir > 0) {
	// erik: what about arrays?
      BEGIN_REFLEVEL_LOOP(cgh) {
	BEGIN_MGLEVEL_LOOP(cgh) {
	  FlipTimeLevels(cgh);
	  for (int rl=0; rl<hh->reflevels(); ++rl) {
	    tt->set_time (rl, mglevel, -1*tt->get_time (rl, mglevel) );
	    tt->set_delta (rl, mglevel,-1*tt->get_delta (rl, mglevel));
	  }
	} END_MGLEVEL_LOOP(cgh);
      } END_REFLEVEL_LOOP(cgh);
	time_dir *= -1;
	cgh->cctk_delta_time = time_dir * base_delta_time / reflevelfact * mglevelfact;
      }
      
      // Evolve each level "backwards" one more timestep
      // Starting with the finest level and proceeding to the coarsest
      BEGIN_REVERSE_REFLEVEL_LOOP(cgh) {
	BEGIN_MGLEVEL_LOOP(cgh) {

          // Restrict from reflevel+1 to reflevel
          // (it checks and does nothing if reflevel is the finest level
          cout << "Calling Restrict for revlevel = " << reflevel
	       << " at time = " << tt->get_time (reflevel, mglevel) << endl;
	  Restrict (cgh);
	  
	  // Cycle time levels
	  // erik: advance time here
	  CycleTimeLevels (cgh);
	  // Advance level times
	  tt->advance_time (reflevel, mglevel);
	  // SHH: What's with this "2*" stuff?
	  Waypoint ("%*sScheduling PRESTEP", 2*reflevel, "");
	  CCTK_ScheduleTraverse ("CCTK_PRESTEP", cgh, CallFunction);
	  Waypoint ("%*sScheduling EVOL", 2*reflevel, "");
	  CCTK_ScheduleTraverse ("CCTK_EVOL", cgh, CallFunction);
	  Waypoint ("%*sScheduling POSTSTEP", 2*reflevel, "");
	  CCTK_ScheduleTraverse ("CCTK_POSTSTEP", cgh, CallFunction);
	  
	  
	} END_MGLEVEL_LOOP(cgh);
      } END_REVERSE_REFLEVEL_LOOP(cgh);
      
      // One final flip to get everything pointed "forward" again
      assert (time_dir < 0);
      // erik: what about arrays?
      BEGIN_MGLEVEL_LOOP(cgh) {
         FlipTimeLevels(cgh);
	 for (int rl=0; rl<hh->reflevels(); ++rl) {
	   tt->set_time (rl, mglevel, 0);
	 }
      } END_MGLEVEL_LOOP(cgh);
      time_dir *= -1;
      cgh->cctk_delta_time = time_dir * base_delta_time / reflevelfact * mglevelfact;

      //SHH:check the current time
      // At this point all "current times" should be the same
      BEGIN_REFLEVEL_LOOP(cgh) {
        BEGIN_MGLEVEL_LOOP(cgh) {
	  cout << "Time at rl=" << reflevel << ", ml=" << mglevel << ", is "
	       << tt->get_time(reflevel, mglevel) << endl;
        } END_MGLEVEL_LOOP(cgh);
      } END_REFLEVEL_LOOP(cgh);

      CCTK_INFO("Finished initializing three timelevels");
    } // end of init_3_timelevels
    
    
    
    BEGIN_REVERSE_REFLEVEL_LOOP(cgh) {
      
      BEGIN_MGLEVEL_LOOP(cgh) {
	
	Waypoint ("%*sCurrent time is %g", 2*reflevel, "", cgh->cctk_time);
	
	// Restrict
	Restrict (cgh);
	
	// Poststep
	Waypoint ("%*sScheduling POSTINITIAL", 2*reflevel, "");
	CCTK_ScheduleTraverse ("CCTK_POSTINITIAL", cgh, CallFunction);
	Waypoint ("%*sScheduling POSTSTEP", 2*reflevel, "");
	CCTK_ScheduleTraverse ("CCTK_POSTSTEP", cgh, CallFunction);
	
	// Checking
	PoisonCheck (cgh, alltimes);
	
	// Checking
	CalculateChecksums (cgh, allbutcurrenttime);
	
	// Recover
	Waypoint ("%*sScheduling RECOVER_VARIABLES", 2*reflevel, "");
	CCTK_ScheduleTraverse ("CCTK_RECOVER_VARIABLES", cgh, CallFunction);
	Waypoint ("%*sScheduling CPINITIAL", 2*reflevel, "");
	CCTK_ScheduleTraverse ("CCTK_CPINITIAL", cgh, CallFunction);
	
	// Analysis
	Waypoint ("%*sScheduling ANALYSIS", 2*reflevel, "");
	CCTK_ScheduleTraverse ("CCTK_ANALYSIS", cgh, CallFunction);
	
	// Output
	Waypoint ("%*sOutputGH", 2*reflevel, "");
	CCTK_OutputGH (cgh);
      
	// Checking
	CheckChecksums (cgh, allbutcurrenttime);
	
      } END_MGLEVEL_LOOP(cgh);
      
    } END_REVERSE_REFLEVEL_LOOP(cgh);
    
    Waypoint ("done with Initialise.");
    
    return 0;
  }
  

} // namespace Carpet
