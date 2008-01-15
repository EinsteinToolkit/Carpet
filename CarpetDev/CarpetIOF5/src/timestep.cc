#include <cassert>
#include <iomanip>
#include <limits>
#include <string>
#include <sstream>

#include <hdf5.h>

#include "cctk.h"
#include "cctk_Functions.h"

#include "timestep.hh"
#include "utils.hh"



namespace CarpetIOF5 {
  
  namespace F5 {
    
    using std::numeric_limits;
    using std::setprecision;
    using std::ostringstream;
    
    
    
    timestep_t::
    timestep_t (file_t & file,
                CCTK_REAL const time,
                char const * const name)
      : m_file (file),
        m_time (time)
    {
      // Construct a group name from the simulation time
      ostringstream buf;
      if (name != 0)
      {
        // Use the name argument if it is present
        buf << name;
      }
      else
      {
        // Create a string from the time without losing information
        int const precision = numeric_limits<CCTK_REAL>::digits10 + 2;
        buf << setprecision (precision) << "t=" << time;
      }
      m_name = buf.str();
      
      m_hdf5_timestep
        = open_or_create_group (m_file.get_hdf5_file(), m_name.c_str());
      assert (m_hdf5_timestep >= 0);
      write_or_check_attribute (m_hdf5_timestep, "time", time);
      
      if (CCTK_IsFunctionAliased ("UniqueSimulationID")) {
        cGH const * const cctkGH = get_file().get_cctkGH();
        char const * const job_id
          = static_cast<char const *> (UniqueSimulationID (cctkGH));
        write_or_check_attribute (m_hdf5_timestep, "simulation id", job_id);
      }
      
      assert (invariant());
    }
    
    
    
    timestep_t::
    ~ timestep_t()
    {
      herr_t const herr = H5Gclose (m_hdf5_timestep);
      assert (not herr);
    }
    
    
      
    file_t & timestep_t::
    get_file ()
      const
    {
      return m_file;
    }
    
    
    
    CCTK_REAL timestep_t::
    get_time ()
      const
    {
      return m_time;
    }
    
    
    
    hid_t timestep_t::
    get_hdf5_timestep()
      const
    {
      return m_hdf5_timestep;
    }
    
    
    
    void timestep_t::
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
        get_file().get_link_destination (l_filename, l_objectname);
        l_objectname += string ("/") + m_name;
      }
      filename = l_filename;
      objectname = l_objectname;
    }
    
    
    
    bool timestep_t::
    invariant()
      const
    {
      return m_hdf5_timestep >= 0;
    }
    
  } // namespace F5

} // namespace CarpetIOF5
