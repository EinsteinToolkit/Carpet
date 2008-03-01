#include <cassert>
#include <cstdio>
#include <cstdlib>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "dist.hh"

#include "carpet.hh"
#include "Timers.hh"



namespace Carpet {
  
  using namespace std;
  
  
  
  int Shutdown (tFleshConfig* fc)
  {
    DECLARE_CCTK_PARAMETERS;
    
    Waypoint ("Starting shutdown");
    
    const int convlev = 0;
    cGH* cctkGH = fc->GH[convlev];
    
    static Timer timer ("Shutdown");
    timer.start();
    for (int rl=reflevels-1; rl>=0; --rl) {
      BEGIN_REVERSE_MGLEVEL_LOOP(cctkGH) {
        enter_level_mode (cctkGH, rl);
        do_global_mode = reflevel==0;
        do_meta_mode = do_global_mode and mglevel==mglevels-1;
        
        Checkpoint ("Shutdown at iteration %d time %g%s%s",
                    cctkGH->cctk_iteration, (double)cctkGH->cctk_time,
                    (do_global_mode ? " (global)" : ""),
                    (do_meta_mode ? " (meta)" : ""));
        
        // Terminate
        Checkpoint ("Scheduling TERMINATE");
        CCTK_ScheduleTraverse ("CCTK_TERMINATE", cctkGH, CallFunction);
        
        leave_level_mode (cctkGH);
      } END_REVERSE_MGLEVEL_LOOP;
    } // for rl
    
    BEGIN_REVERSE_MGLEVEL_LOOP(cctkGH) {
      do_global_mode = true;
      do_meta_mode = mglevel==mglevels-1;
      
      // Shutdown
      Checkpoint ("Scheduling SHUTDOWN");
      CCTK_ScheduleTraverse ("CCTK_SHUTDOWN", cctkGH, CallFunction);
      
    } END_REVERSE_MGLEVEL_LOOP;
    timer.stop();
    if (output_timers_every > 0) {
      TimerSet::writeData (cctkGH, timer_file);
    }
    
    // earlier checkpoint before finalising MPI
    Waypoint ("Done with shutdown");
    
    return 0;
  }
  
} // namespace Carpet
