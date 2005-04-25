#ifndef EXTENDING_HH
#define EXTENDING_HH

#include <set>

#include "cctk.h"



namespace CarpetIOF5 {
  
  using std::set;
  
  
  
  class extending {
    
    static char const * const extension_name;
    
    struct extension_t {
      set<char const *> did_truncate;
    };
    
    extension_t * m_extension;
    
  public:
    
    static void
    create (cGH * cctkGH);
    
  private:
    
    static void *
    setup (tFleshConfig * config,
           int convlevel,
           cGH * cctkGH);
    
  public:
    
    extending (cGH const * cctkGH);
    
    set<char const *>
    get_did_truncate ()
      const;
    
  };
  
} // namespace CarpetIOF5

#endif  // #ifndef EXDENDING_HH
