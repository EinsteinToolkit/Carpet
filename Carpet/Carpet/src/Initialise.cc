#include <assert.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctki_GHExtensions.h"
#include "cctki_ScheduleBindings.h"
#include "cctki_WarnLevel.h"

#include "carpet.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Initialise.cc,v 1.11 2002/06/06 19:48:29 shawley Exp $";

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
    Waypoint ("PARAMCHECK");
    CCTK_ScheduleTraverse ("CCTK_PARAMCHECK", cgh, CallFunction);
    CCTKi_FinaliseParamWarn();
    set_mglevel (cgh, -1);
    
    Waypoint ("Initialising iteration %d...", cgh->cctk_iteration);
    
    BEGIN_REFLEVEL_LOOP(cgh) {
      
      BEGIN_MGLEVEL_LOOP(cgh) {
	
	// Checking
	Poison (cgh, alltimes);
	
	// Set up the grid
	Waypoint ("%*sScheduling BASEGRID", 2*reflevel, "");
	CCTK_ScheduleTraverse ("CCTK_BASEGRID", cgh, CallFunction);
	if (reflevel==0) {
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


    // Scott's algorithm for getting two extra timelevels of data
    if (init_3_timelevels) {
       int time_dir = 1; //Positive = forward (+t), Negative = backward (-t)
       CCTK_INFO("Initializing THREE timelevels");
       
       BEGIN_REFLEVEL_LOOP(cgh) {
          BEGIN_MGLEVEL_LOOP(cgh) {
	      // Evolve "forward" (which may be backward for lev=1,3,5,7...)
	      // Cycle time levels
	      CycleTimeLevels (cgh);

	      Waypoint ("%*sScheduling PRESTEP", 2*reflevel, "");
	      CCTK_ScheduleTraverse ("CCTK_PRESTEP", cgh, CallFunction);
	      Waypoint ("%*sScheduling EVOL", 2*reflevel, "");
	      CCTK_ScheduleTraverse ("CCTK_EVOL", cgh, CallFunction);
	      Waypoint ("%*sScheduling POSTSTEP", 2*reflevel, "");
	      CCTK_ScheduleTraverse ("CCTK_POSTSTEP", cgh, CallFunction);

	      FlipTimeLevels(cgh);
              cgh->cctk_delta_time *= -1;
              //    Keep track of which direction (in time) we're integrating
	      time_dir *= -1;
	
	      // Evolve in the opposite time-direction
	      // Cycle time levels
	      CycleTimeLevels (cgh);
	      Waypoint ("%*sScheduling PRESTEP", 2*reflevel, "");
	      CCTK_ScheduleTraverse ("CCTK_PRESTEP", cgh, CallFunction);
	      Waypoint ("%*sScheduling EVOL", 2*reflevel, "");
	      CCTK_ScheduleTraverse ("CCTK_EVOL", cgh, CallFunction);
	      Waypoint ("%*sScheduling POSTSTEP", 2*reflevel, "");
	      CCTK_ScheduleTraverse ("CCTK_POSTSTEP", cgh, CallFunction);

          } END_MGLEVEL_LOOP(cgh);
       } END_REFLEVEL_LOOP(cgh);


       // Make sure we're pointed backwards, in order to get 2 "previous"
       //   timelevels.  We could change the if statement to 
       //   "if (mod(MaxLevels,2) == 0)", but I prefer to check time_dir 
       //   explicitly, because it's easier to follow and I don't have to
       //   worry about having made a mistake
       if (time_dir > 0) {
	      FlipTimeLevels(cgh);
              cgh->cctk_delta_time *= -1;
	      time_dir *= -1;
       }

       // Evolve each level backwards one more timestep
       BEGIN_REVERSE_REFLEVEL_LOOP(cgh) {
           BEGIN_MGLEVEL_LOOP(cgh) {

	      // Cycle time levels
	      CycleTimeLevels (cgh);
	      Waypoint ("%*sScheduling PRESTEP", 2*reflevel, "");
	      CCTK_ScheduleTraverse ("CCTK_PRESTEP", cgh, CallFunction);
	      Waypoint ("%*sScheduling EVOL", 2*reflevel, "");
	      CCTK_ScheduleTraverse ("CCTK_EVOL", cgh, CallFunction);
	      Waypoint ("%*sScheduling POSTSTEP", 2*reflevel, "");
	      CCTK_ScheduleTraverse ("CCTK_POSTSTEP", cgh, CallFunction);

	      // Restrict
	      Restrict (cgh);


           } END_MGLEVEL_LOOP(cgh);
       } END_REVERSE_REFLEVEL_LOOP(cgh);

       CCTK_INFO("Finished initializing three timelevels");
    }  // end of init_3_timelevels

    
    BEGIN_REVERSE_REFLEVEL_LOOP(cgh) {
      
      BEGIN_MGLEVEL_LOOP(cgh) {
	
	// Restrict
	Restrict (cgh);
	
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
