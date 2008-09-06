#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <sstream>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "dist.hh"

#include "carpet.hh"
#include "Timers.hh"



namespace Carpet {
  
  using namespace std;
  
  
  
  int OutputGH (cGH const * const cctkGH)
  {
    DECLARE_CCTK_PARAMETERS;
    
    static Timer timer ("OutputGH");
    timer.start();
    
    Checkpoint ("OutputGH");
    
    int const num_methods = CCTK_NumIOMethods();
    if (num_methods == 0) {
      return -1;
    }
    
    static vector<Timer *> timers;
    timers.resize (num_methods, NULL);
    
    int num_vars = 0;
    for (int handle = 0; handle < num_methods; ++ handle) {
      
      IOMethod const * const method = CCTK_IOMethod (handle);
      assert (method);
      
      if (not timers.at(handle)) {
        ostringstream buf;
        buf << "OutputGH"
            << "::" << method->implementation
            << "::" << method->name
            << " [" << handle << "]";
        timers.at(handle) = new Timer (buf.str().c_str());
      }
      
      timers.at(handle)->start();
      num_vars += method->OutputGH (cctkGH);
      timers.at(handle)->stop();
      
    } // for handle
    
    timer.stop();
    
    return num_vars;
  }
  
} // namespace Carpet
