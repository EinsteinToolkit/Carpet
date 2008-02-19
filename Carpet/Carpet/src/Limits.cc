#include <cctk.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <iostream>
#include <sys/resource.h>

namespace Carpet {
  
  using namespace std;
  
  
  
  static
  void
  set_limit (int resource, char const * name, CCTK_INT value);
  
  static
  ostream &
  operator<< (ostream & s, struct rlimit const & limit);
  
  static
  void
  output (ostream & s, rlim_t const & value);
  
  
  
  void
  SetSystemLimits ()
  {
    DECLARE_CCTK_PARAMETERS;
    set_limit (RLIMIT_CORE, "core file size", max_core_size_MB);
    set_limit (RLIMIT_AS,   "memory size",    max_memory_size_MB);
  }
  
  
  
  void
  set_limit (int const resource, char const * const name, CCTK_INT const value)
  {
    struct rlimit limit;
    int const ierr1 = getrlimit (resource, & limit);
    assert (not ierr1);
    
    if (value == -2 ) {
      // Only show limit
      cout << "Current " << name << " limit: " << limit << endl;
      return;
    }
    
    cout << "Old " << name << " limit: " << limit << endl;
    
    if (value == -1) {
      limit.rlim_cur = limit.rlim_max;
    } else {
      limit.rlim_cur = min ((rlim_t) value * 1024 * 1024, limit.rlim_max);
    }
    
    int const ierr2 = setrlimit (resource, & limit);
    assert (not ierr2);
    
    int const ierr3 = getrlimit (resource, & limit);
    assert (not ierr3);
    
    cout << "New " << name << " limit: " << limit << endl;
  }
  
  
  
  static
  ostream &
  operator<< (ostream & s, struct rlimit const & limit)
  {
    s << "hard=";
    output (s, limit.rlim_max);
    s << ", soft=";
    output (s, limit.rlim_cur);
    return s;
  }
  
  
  
  static
  void
  output (ostream & s, rlim_t const & value)
  {
    if (value == RLIM_INFINITY) {
      s << "[unlimited]";
    } else {
      s << ((CCTK_REAL) value / 1.0e+6) << " MB";
    }
  }
  
} // namespace Carpet
