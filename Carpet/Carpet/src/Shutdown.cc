#include <cassert>
#include <cstdio>
#include <cstdlib>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "dist.hh"

#include "carpet.hh"



namespace Carpet {
  
  using namespace std;
  
  
  
  int Shutdown (tFleshConfig* fc)
  {
    DECLARE_CCTK_PARAMETERS;
    
    Waypoint ("Starting shutdown");
    
    const int convlev = 0;
    cGH* cgh = fc->GH[convlev];
    
    for (int rl=reflevels-1; rl>=0; --rl) {
      BEGIN_REVERSE_MGLEVEL_LOOP(cgh) {
        enter_level_mode (cgh, rl);
        do_global_mode = reflevel==0;
        do_meta_mode = do_global_mode and mglevel==mglevels-1;
        
        Checkpoint ("Shutdown at iteration %d time %g%s%s",
                    cgh->cctk_iteration, (double)cgh->cctk_time,
                    (do_global_mode ? " (global)" : ""),
                    (do_meta_mode ? " (meta)" : ""));
        
        // Terminate
        Checkpoint ("Scheduling TERMINATE");
        CCTK_ScheduleTraverse ("CCTK_TERMINATE", cgh, CallFunction);
        
        leave_level_mode (cgh);
      } END_REVERSE_MGLEVEL_LOOP;
    } // for rl
    
    BEGIN_REVERSE_MGLEVEL_LOOP(cgh) {
      do_global_mode = true;
      do_meta_mode = mglevel==mglevels-1;
      
      // Shutdown
      Checkpoint ("Scheduling SHUTDOWN");
      CCTK_ScheduleTraverse ("CCTK_SHUTDOWN", cgh, CallFunction);
      
    } END_REVERSE_MGLEVEL_LOOP;
    
    CCTK_PRINTSEPARATOR;
    printf ("Done.\n");
    
    // earlier checkpoint before finalising MPI
    Waypoint ("Done with shutdown");
    
    dist::finalize();
    
    return 0;
  }
  
} // namespace Carpet
