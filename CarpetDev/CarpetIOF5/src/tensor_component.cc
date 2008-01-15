#include <cassert>
#include <cstdlib>

#include "cctk.h"

#include "tensor_component.hh"



namespace CarpetIOF5 {
  
  namespace F5 {
    
    tensor_component_t::
    tensor_component_t (physical_quantity_t & physical_quantity,
                        int const variable)
      : m_physical_quantity (physical_quantity),
        m_variable (variable)
    {
      assert (variable >= 0 and variable < CCTK_NumVars());
      
      char const * const name = CCTK_VarName (variable);
      assert (name != 0);
      m_name = string (name);
      
      m_hdf5_tensor_component
        = open_or_create_group (m_physical_quantity
                                .get_hdf5_physical_quantity(),
                                name);
      assert (m_hdf5_tensor_component >= 0);
      
      assert (invariant());
    }
    
    
    
    tensor_component_t::
    ~ tensor_component_t ()
    {
      herr_t const herr = H5Gclose (m_hdf5_tensor_component);
      assert (not herr);
    }
    
    
      
    physical_quantity_t & tensor_component_t::
    get_physical_quantity ()
        const
    {
      return m_physical_quantity;
    }
    
    
    
    hid_t tensor_component_t::
    get_variable ()
      const
    {
      return m_variable;
    }
    
    
    
    hid_t tensor_component_t::
    get_hdf5_tensor_component ()
      const
    {
      return m_hdf5_tensor_component;
    }
    
    
    
    void tensor_component_t::
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
        get_physical_quantity().get_link_destination (l_filename, l_objectname);
        l_objectname += string ("/") + m_name;
      }
      filename = l_filename;
      objectname = l_objectname;
    }
    
    
    
    bool tensor_component_t::
    invariant ()
      const
    {
      return (m_variable >= 0 and m_variable < CCTK_NumVars()
              and m_hdf5_tensor_component >= 0);
    }
    
  } // namespace F5

} // namespace CarpetIOF5
