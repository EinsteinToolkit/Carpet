#include <cassert>
#include <cstdlib>
#include <sstream>
#include <string>

#include "cctk.h"

#include "carpet.hh"

#include "data_region.hh"
#include "meta_data_region.hh"
#include "utils.hh"



namespace CarpetIOF5 {
  
  namespace F5 {
    
    using std::ostringstream;
    
    
    
    meta_data_region_t::
    meta_data_region_t (tensor_component_t & tensor_component,
                        bbox<int, dim> const & region)
      : m_tensor_component (tensor_component),
        m_region (region)
    {
      assert (not region.empty());
      
      string const namestr = data_region_t::name_from_region (region);
      char const * const name = namestr.c_str();
      assert (name != 0);
      
      assert (invariant());
    }
    
    
    
    meta_data_region_t::
    ~ meta_data_region_t ()
    {
    }
    
    
    
    tensor_component_t & meta_data_region_t::
    get_tensor_component ()
      const
    {
      return m_tensor_component;
    }
    
    
    
    void meta_data_region_t::
    write (int const proc)
      const
    {
      string filename;
      string objectname;
      get_tensor_component().get_link_destination (filename, objectname);
      
      string const name = data_region_t::name_from_region (m_region);
      
      hid_t const hdf5_tensor_component
        = m_tensor_component.get_hdf5_tensor_component();
      
      herr_t const herr
        = H5Lcreate_external (filename.c_str(),
                              objectname.c_str(),
                              hdf5_tensor_component,
                              name.c_str(),
                              H5P_DEFAULT, H5P_DEFAULT);
      assert (not herr);
    }
    
    
    
    bool meta_data_region_t::
    invariant ()
      const
    {
      return not m_region.empty();
    }
    
  } // namespace F5

} // namespace CarpetIOF5
