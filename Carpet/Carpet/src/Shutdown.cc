#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "dist.hh"

#include "carpet.hh"

extern "C" {
  static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Shutdown.cc,v 1.12 2003/08/03 17:09:02 schnetter Exp $";
  CCTK_FILEVERSION(Carpet_Carpet_Shutdown_cc);
}



namespace Carpet {
  
  using namespace std;
  
  
  
  int Shutdown (tFleshConfig* fc)
  {
    DECLARE_CCTK_PARAMETERS;
    
    Waypoint ("starting Shutdown...");
    
    const int convlev = 0;
    cGH* cgh = fc->GH[convlev];
    
    Waypoint ("Current time is %g", cgh->cctk_time);
    
    BEGIN_REFLEVEL_LOOP(cgh) {
      BEGIN_MGLEVEL_LOOP(cgh) {
	
        do_global_mode = reflevel == 0;
        
	Waypoint ("%*sCurrent time is %g%s", 2*reflevel, "",
                  cgh->cctk_time,
                  do_global_mode ? "   (global time)" : "");
        
        // Terminate
        Waypoint ("%*sScheduling TERMINATE", 2*reflevel, "");
        CCTK_ScheduleTraverse ("CCTK_TERMINATE", cgh, CallFunction);
	
      } END_MGLEVEL_LOOP;
    } END_REFLEVEL_LOOP;
    
    do_global_mode = true;
    
    // Shutdown
    Waypoint ("Scheduling SHUTDOWN");
    CCTK_ScheduleTraverse ("CCTK_SHUTDOWN", cgh, CallFunction);
    
    CCTK_PRINTSEPARATOR;
    printf ("Done.\n");
    
    // earlier checkpoint before finalising MPI
    Waypoint ("done with Shutdown.");
    
    dist::finalize();
    
    return 0;
  }
  
} // namespace Carpet
