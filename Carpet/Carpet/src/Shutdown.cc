#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "Carpet/CarpetLib/src/dist.hh"

#include "carpet.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Shutdown.cc,v 1.6 2002/03/26 13:22:27 schnetter Exp $";

CCTK_FILEVERSION(Carpet_Shutdown_cc)



namespace Carpet {
  
  using namespace std;
  
  
  
  int Shutdown (tFleshConfig* fc)
  {
    DECLARE_CCTK_PARAMETERS;
    
    Waypoint ("starting Shutdown...");
    
    const int convlev = 0;
    cGH* cgh = fc->GH[convlev];
    
    set_mglevel (cgh, 0);

    // Terminate
    Waypoint ("%*sScheduling TERMINATE", 2*reflevel, "");
    CCTK_ScheduleTraverse ("CCTK_TERMINATE", cgh, CallFunction);
    
    // Shutdown
    Waypoint ("%*sScheduling SHUTDOWN", 2*reflevel, "");
    CCTK_ScheduleTraverse ("CCTK_SHUTDOWN", cgh, CallFunction);
    
    set_mglevel (cgh, -1);
    
    CCTK_PRINTSEPARATOR;
    printf ("Done.\n");
    
    // earlier checkpoint before calling finalising MPI
    Waypoint ("done with Shutdown.");
    
    dist::finalize();
    
    return 0;
  }
  
} // namespace Carpet
