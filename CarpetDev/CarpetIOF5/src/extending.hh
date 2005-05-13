#ifndef EXTENDING_HH
#define EXTENDING_HH

#include <set>
#include <vector>

#include "cctk.h"



namespace CarpetIOF5 {
  
  using std::set;
  using std::vector;
  
  
  
  class extending_t {
    
    static char const * const extension_name;
    
    struct extension_t {
      set<char const *> did_truncate;
      // [mglevel][reflevel][variable];
      vector<vector<vector<int> > > last_output_iteration;
      vector<vector<vector<CCTK_REAL> > > last_output_time;
    };
    
    extension_t * m_extension;
    
  public:
    
    static void
    create ();
    
    static void *
    setup (tFleshConfig * config,
           int convlevel,
           cGH * cctkGH);
    
    extending_t (cGH const * cctkGH);
    
    bool
    get_did_truncate (char const * name)
      const;
    
    void
    set_did_truncate (char const * name);
    
    int
    get_last_output_iteration (int ml, int rl, int vi)
      const;
    
    void
    set_last_output_iteration (int ml, int rl, int vi, int iteration);
    
    CCTK_REAL
    get_last_output_time (int ml, int rl, int vi)
      const;
    
    void
    set_last_output_time (int ml, int rl, int vi, CCTK_REAL time);
    
  private:
    
    template<typename T>
    static void
    resize_last_output (int ml, int rl, int vi,
                        vector<vector<vector<T> > > & array);
    
  };
  
} // namespace CarpetIOF5

#endif  // #ifndef EXDENDING_HH
