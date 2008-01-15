#ifndef META_DATA_REGION_HH
#define META_DATA_REGION_HH

#include <hdf5.h>

#include "bbox.hh"
#include "defs.hh"

#include "tensor_component.hh"



namespace CarpetIOF5 {
  
  namespace F5 {
    
    class meta_data_region_t {
      
      tensor_component_t & m_tensor_component;
      
      bbox<int, dim> const m_region;
      
      hid_t m_properties;
      hid_t m_dataset;
      hid_t m_dataspace;
      
      meta_data_region_t ();
      meta_data_region_t (meta_data_region_t const &);
      meta_data_region_t operator= (meta_data_region_t const &);
      
    public:
      
      meta_data_region_t (tensor_component_t & tensor_component,
                          bbox<int, dim> const & region);
      
      virtual
      ~ meta_data_region_t ();
      
      tensor_component_t &
      get_tensor_component ()
        const;
      
      void
      write (int proc)
        const;
      
      virtual bool
      invariant ()
        const;
      
    };
    
  } // namespace F5

} // namespace CarpetIOF5

#endif  // #ifndef META_DATA_REGION_HH
