#include <cassert>

#include <hdf5.h>

#include "cctk.h"

#include "simulation.hh"
#include "utils.hh"



namespace CarpetIOF5 {
  
  namespace F5 {
    
    using std::ostringstream;
    
    
    
    simulation_t::
    simulation_t (timestep_t & timestep,
                  char const * const name)
      : m_timestep (timestep),
        m_name (name)
    {
      m_hdf5_simulation
        = open_or_create_group (m_timestep.get_hdf5_timestep(),
                                m_name.c_str());
      assert (m_hdf5_simulation >= 0);
      
      assert (invariant());
    }
    
    
    
    simulation_t::
    ~ simulation_t()
    {
      herr_t const herr = H5Gclose (m_hdf5_simulation);
      assert (not herr);
    }
    
    
    
    timestep_t & simulation_t::
    get_timestep ()
      const
    {
      return m_timestep;
    }
    
    
    
    hid_t simulation_t::
    get_hdf5_simulation()
      const
    {
      return m_hdf5_simulation;
    }
    
    
    
    void simulation_t::
    get_link_destination (string & filename,
                          string & objectname)
      const
    {
      static bool initialised = false;
      static string l_filename;
      static string l_objectname;
      if (not initialised)
      {
        initialised = true;
        get_timestep().get_link_destination (l_filename, l_objectname);
        l_objectname += string ("/") + m_name;
      }
      filename = l_filename;
      objectname = l_objectname;
    }
    
    
    
    bool simulation_t::
    invariant()
      const
    {
      return m_hdf5_simulation >= 0;
    }
    
  } // namespace F5

} // namespace CarpetIOF5
