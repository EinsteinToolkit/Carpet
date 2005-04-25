#ifndef SIMULATION_HH
#define SIMULATION_HH

#include <set>
#include <string>

#include <hdf5.h>

#include "cctk.h"

#include "timestep.hh"



namespace CarpetIOF5 {
  
  using std::string;
  
  namespace F5 {
  
    class simulation_t {
      
      timestep_t & m_timestep;
      
      string m_name;
      
      hid_t m_hdf5_simulation;
      
    public:
      
      simulation_t (timestep_t & timestep,
                    char const * name);
      
      virtual
      ~ simulation_t ();
      
      hid_t
      get_hdf5_simulation ()
        const;
      
      virtual bool
      invariant ()
        const;
    };
    
  } // namespace F5

} // namespace CarpetIOF5

#endif  // #ifndef SIMULATION_HH
