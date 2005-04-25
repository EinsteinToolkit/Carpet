#ifndef TIMESTEP_HH
#define TIMESTEP_HH

#include <string>

#include <hdf5.h>

#include "cctk.h"

#include "file.hh"



namespace CarpetIOF5 {
  
  using std::string;
  
  namespace F5 {
  
    class timestep_t {
      
      file_t & m_file;
      
      CCTK_REAL const m_time;
      
      string m_name;
      
      hid_t m_hdf5_timestep;
      
    public:
      
      timestep_t (file_t & file,
                  CCTK_REAL time,
                  char const * name = 0);
      
      virtual
      ~ timestep_t ();
      
      hid_t
      get_hdf5_timestep ()
        const;
      
      virtual bool
      invariant ()
        const;
    };
    
  } // namespace F5

} // namespace CarpetIOF5

#endif  // #ifndef TIMESTEP_HH
