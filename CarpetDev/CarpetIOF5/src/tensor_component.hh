#ifndef TENSOR_COMPONENT_HH
#define TENSOR_COMPONENT_HH

#include <hdf5.h>

#include "physical_quantity.hh"



namespace CarpetIOF5 {
  
  namespace F5 {
    
    class tensor_component_t {
      
      physical_quantity_t & m_physical_quantity;
      
      int const m_variable;
      
      hid_t m_hdf5_tensor_component;
      
    public:
      
      tensor_component_t (physical_quantity_t & physical_quantity,
                          int variable);
      
      virtual
      ~ tensor_component_t ();
      
      hid_t
      get_variable ()
        const;
      
      hid_t
      get_hdf5_tensor_component ()
        const;
      
      virtual bool
      invariant ()
        const;
      
    };
    
  } // namespace F5

} // namespace CarpetIOF5

#endif  // #ifndef TENSOR_COMPONENT_HH
