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
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Initialise.cc,v 1.24 2002/12/12 12:55:45 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_Carpet_Initialise_cc);
}



namespace Carpet {
  
  using namespace std;
  
  
  
  void Initialise3TL (cGH* const cgh);
  
  
  
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
	Waypoint ("%*sScheduling POSTINITIAL", 2*reflevel, "");
	CCTK_ScheduleTraverse ("CCTK_POSTINITIAL", cgh, CallFunction);
	
	// Poststep
	Waypoint ("%*sScheduling POSTSTEP", 2*reflevel, "");
	CCTK_ScheduleTraverse ("CCTK_POSTSTEP", cgh, CallFunction);
	
	// Checking
	PoisonCheck (cgh, alltimes);
	
      } END_MGLEVEL_LOOP(cgh);
      
      // Regrid
      Waypoint ("%*sRegrid", 2*reflevel, "");
      Regrid (cgh);
      
    } END_REFLEVEL_LOOP(cgh);
    
    
    
    BEGIN_REVERSE_REFLEVEL_LOOP(cgh) {
      
      BEGIN_MGLEVEL_LOOP(cgh) {
	
	// Restrict
	Waypoint ("%*sCurrent time is %g", 2*reflevel, "", cgh->cctk_time);
	Restrict (cgh);
	
      } END_MGLEVEL_LOOP(cgh);
      
    } END_REVERSE_REFLEVEL_LOOP(cgh);
    
    
    
    if (init_3_timelevels) {
      // Call Scott's algorithm for getting two extra timelevels of data
      Initialise3TL (cgh);
    }
    
    
    
    BEGIN_REFLEVEL_LOOP(cgh) {
      
      BEGIN_MGLEVEL_LOOP(cgh) {
	
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
      
    } END_REFLEVEL_LOOP(cgh);
    
    Waypoint ("done with Initialise.");
    
    return 0;
  }
  
  
  
  void Initialise3TL (cGH* const cgh)
  {
    DECLARE_CCTK_PARAMETERS;
    
    CCTK_INFO ("Initialising three timelevels");
    
    // Act as if this was the first iteration
    cgh->cctk_iteration = 1;
    
    BEGIN_REFLEVEL_LOOP(cgh) {
      BEGIN_MGLEVEL_LOOP(cgh) {
	
	cgh->cctk_time = cctk_initial_time;
	
	// Cycle time levels (ignore arrays)
	cout << "3TL rl=" << reflevel << " cycling" << endl;
	CycleTimeLevels (cgh);
	
	// Advance level times
	tt->advance_time (reflevel, mglevel);
	cgh->cctk_time += cgh->cctk_delta_time;
	cout << "3TL rl=" << reflevel << " ml=" << mglevel
	     << " time=" << tt->get_time (reflevel, mglevel)
	     << " time=" << cgh->cctk_time / base_delta_time << endl;
	
	// Evolve forward
	Waypoint ("%*sScheduling PRESTEP", 2*reflevel, "");
	CCTK_ScheduleTraverse ("CCTK_PRESTEP", cgh, CallFunction);
	Waypoint ("%*sScheduling EVOL", 2*reflevel, "");
	CCTK_ScheduleTraverse ("CCTK_EVOL", cgh, CallFunction);
	Waypoint ("%*sScheduling POSTSTEP", 2*reflevel, "");
	CCTK_ScheduleTraverse ("CCTK_POSTSTEP", cgh, CallFunction);
	
	// Flip time levels
	cout << "3TL rl=" << reflevel << " flipping" << endl;
	FlipTimeLevels (cgh);
	
	// Invert level times
	for (int rl=0; rl<hh->reflevels(); ++rl) {
	  tt->set_delta (rl, mglevel, - tt->get_delta (rl, mglevel));
	  tt->advance_time (rl, mglevel);
	  tt->advance_time (rl, mglevel);
	}
	cgh->cctk_delta_time *= -1;
	cgh->cctk_time += cgh->cctk_delta_time;
	cgh->cctk_time += cgh->cctk_delta_time;
	cout << "3TL rl=" << reflevel << " ml=" << mglevel
	     << " time=" << tt->get_time (reflevel, mglevel)
	     << " time=" << cgh->cctk_time / base_delta_time << endl;
	
	// Evolve backward
	Waypoint ("%*sScheduling PRESTEP", 2*reflevel, "");
	CCTK_ScheduleTraverse ("CCTK_PRESTEP", cgh, CallFunction);
	Waypoint ("%*sScheduling EVOL", 2*reflevel, "");
	CCTK_ScheduleTraverse ("CCTK_EVOL", cgh, CallFunction);
	Waypoint ("%*sScheduling POSTSTEP", 2*reflevel, "");
	CCTK_ScheduleTraverse ("CCTK_POSTSTEP", cgh, CallFunction);
	
	// Flip time levels back
	cout << "3TL rl=" << reflevel << " flipping back" << endl;
	FlipTimeLevels (cgh);
	
	// Invert level times back
	for (int rl=0; rl<hh->reflevels(); ++rl) {
	  tt->set_delta (rl, mglevel, - tt->get_delta (rl, mglevel));
	  tt->advance_time (rl, mglevel);
	  tt->advance_time (rl, mglevel);
	}
	cgh->cctk_delta_time *= -1;
	cgh->cctk_time += cgh->cctk_delta_time;
	cgh->cctk_time += cgh->cctk_delta_time;
	cout << "3TL rl=" << reflevel << " ml=" << mglevel
	     << " time=" << tt->get_time (reflevel, mglevel)
	     << " time=" << cgh->cctk_time / base_delta_time << endl;
	
      } END_MGLEVEL_LOOP(cgh);
    } END_REFLEVEL_LOOP(cgh);
    
    cout << "Hourglass structure in place" << endl;
      
    // Act as if this was the second iteration
    cgh->cctk_iteration = 2;
    
    // Evolve each level "backwards" one more timestep
    // Starting with the finest level and proceeding to the coarsest
    BEGIN_REVERSE_REFLEVEL_LOOP(cgh) {
      BEGIN_MGLEVEL_LOOP(cgh) {
	
	cgh->cctk_time = cctk_initial_time + cgh->cctk_delta_time;
	
	// Restrict
	cout << "3TL rl=" << reflevel << " restricting" << endl;
	Restrict (cgh);
	
	// Flip time levels
	cout << "3TL rl=" << reflevel << " flipping" << endl;
	FlipTimeLevels (cgh);
	
	// Invert level times
	for (int rl=0; rl<hh->reflevels(); ++rl) {
	  tt->set_delta (rl, mglevel, - tt->get_delta (rl, mglevel));
	  tt->advance_time (rl, mglevel);
	  tt->advance_time (rl, mglevel);
	}
	cgh->cctk_delta_time *= -1;
	cgh->cctk_time += cgh->cctk_delta_time;
	cgh->cctk_time += cgh->cctk_delta_time;
	cout << "3TL rl=" << reflevel << " ml=" << mglevel
	     << " time=" << tt->get_time (reflevel, mglevel)
	     << " time=" << cgh->cctk_time / base_delta_time << endl;
	
	// Cycle time levels
	cout << "3TL rl=" << reflevel << " cycling" << endl;
	CycleTimeLevels (cgh);
	
	// Advance level times
	tt->advance_time (reflevel, mglevel);
	cgh->cctk_time += cgh->cctk_delta_time;
	cout << "3TL rl=" << reflevel << " ml=" << mglevel
	     << " time=" << tt->get_time (reflevel, mglevel)
	     << " time=" << cgh->cctk_time / base_delta_time << endl;
	
	// Evolve backward
	Waypoint ("%*sScheduling PRESTEP", 2*reflevel, "");
	CCTK_ScheduleTraverse ("CCTK_PRESTEP", cgh, CallFunction);
	Waypoint ("%*sScheduling EVOL", 2*reflevel, "");
	CCTK_ScheduleTraverse ("CCTK_EVOL", cgh, CallFunction);
	Waypoint ("%*sScheduling POSTSTEP", 2*reflevel, "");
	CCTK_ScheduleTraverse ("CCTK_POSTSTEP", cgh, CallFunction);
	
	// Flip time levels back
	cout << "3TL rl=" << reflevel << " flipping back" << endl;
	FlipTimeLevels (cgh);
	
	// Invert level times back
	for (int rl=0; rl<hh->reflevels(); ++rl) {
	  tt->set_delta (rl, mglevel, - tt->get_delta (rl, mglevel));
	  tt->advance_time (rl, mglevel);
	  tt->advance_time (rl, mglevel);
	}
	cgh->cctk_delta_time *= -1;
	cgh->cctk_time += cgh->cctk_delta_time;
	cgh->cctk_time += cgh->cctk_delta_time;
	cout << "3TL rl=" << reflevel << " ml=" << mglevel
	     << " time=" << tt->get_time (reflevel, mglevel)
	     << " time=" << cgh->cctk_time / base_delta_time << endl;
	
      } END_MGLEVEL_LOOP(cgh);
    } END_REVERSE_REFLEVEL_LOOP(cgh);
    
    // Reset stuff
    cgh->cctk_iteration = 0;
    cgh->cctk_time = cctk_initial_time;
    
    CCTK_INFO ("Finished initialising three timelevels");
  }
  
} // namespace Carpet
