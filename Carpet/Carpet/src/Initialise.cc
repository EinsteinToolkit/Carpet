#include <assert.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctki_GHExtensions.h"
#include "cctki_ScheduleBindings.h"
#include "cctki_WarnLevel.h"

#include "carpet.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Initialise.cc,v 1.7 2002/01/14 14:59:24 schnetter Exp $";



namespace Carpet {
  
  using namespace std;
  
  
  
  int Initialise (tFleshConfig* fc)
  {
    cout << "zwei..." << endl;
    DECLARE_CCTK_PARAMETERS;
    cout << "zwei..." << endl;
    
    // Initialise stuff
    const int convlev = 0;
    cout << "zwei..." << endl;
    cGH* const cgh = CCTK_SetupGH (fc, convlev);
    cout << "zwei..." << endl;
    CCTKi_AddGH (fc, convlev, cgh);
    cout << "zwei..." << endl;
    
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
