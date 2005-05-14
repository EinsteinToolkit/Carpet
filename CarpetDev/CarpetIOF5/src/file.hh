#ifndef FILE_HH
#define FILE_HH

#include <string>

#include <hdf5.h>

#include "cctk.h"

#include "extending.hh"



namespace CarpetIOF5 {
  
  using std::string;
  
  
  
  namespace F5 {
    
    class file_t {
      
      cGH const * const m_cctkGH;
      
      string const m_filename;
      
      hid_t m_hdf5_file;
      
      file_t ();
      file_t (file_t const &);
      file_t operator= (file_t const &);
      
    public:
      
      file_t (cGH const * cctkGH,
              string filename,
              bool do_truncate);
      
      virtual
      ~ file_t ();
      
      cGH const *
      get_cctkGH ()
        const;
      
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
