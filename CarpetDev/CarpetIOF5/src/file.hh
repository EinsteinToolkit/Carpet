#ifndef FILE_HH
#define FILE_HH

#include <hdf5.h>

#include "cctk.h"

#include "extending.hh"



namespace CarpetIOF5 {
  
  namespace F5 {
    
    class file_t {
      
      cGH const * const m_cctkGH;
      
      char const * const m_filename;
      
      hid_t m_hdf5_file;
      
    public:
      
      file_t (cGH const * cctkGH,
              char const * filename);
      
      virtual
      ~ file_t ();
      
      hid_t
      get_hdf5_file ()
        const;
      
      virtual bool
      invariant ()
        const;
    };
    
  } // namespace F5

} // namespace CarpetIOF5

#endif  // #ifndef FILE_HH
