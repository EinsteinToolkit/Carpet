#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "Carpet/CarpetLib/src/dist.hh"

#include "carpet.hh"

static const char* rcsid = "$Header: /home/eschnett/C/carpet/Carpet/Carpet/Carpet/src/Shutdown.cc,v 1.1 2001/07/04 12:29:47 schnetter Exp $";



namespace Carpet {
  
  using namespace std;
  
  
  
  int Shutdown (tFleshConfig* fc)
  {
    DECLARE_CCTK_PARAMETERS;
    
    Checkpoint ("starting Shutdown...");
    
    const int convlev = 0;
    cGH* cgh = fc->GH[convlev];
    
    // Terminate
    BEGIN_REFLEVEL_LOOP(cgh) {
      Checkpoint ("%*sScheduling TERMINATE", 2*reflevel, "");
      CCTK_ScheduleTraverse ("CCTK_TERMINATE", cgh, CallFunction);
    } END_REFLEVEL_LOOP(cgh);
    
    // Shutdown
    BEGIN_REFLEVEL_LOOP(cgh) {
      Checkpoint ("%*sScheduling SHUTDOWN", 2*reflevel, "");
      CCTK_ScheduleTraverse ("CCTK_SHUTDOWN", cgh, CallFunction);
    } END_REFLEVEL_LOOP(cgh);
    
    CCTK_PRINTSEPARATOR;
    printf ("Done.\n");
    
    // earlier checkpoint before calling finalising MPI
    Checkpoint ("done with Shutdown.");
    
    dist::finalize();
    
    return 0;
  }
  
} // namespace Carpet
