#ifndef DATA_REGION_HH
#define DATA_REGION_HH

// force HDF5 1.8.x installations to use the new API
#ifdef H5Dcreate_vers
#undef H5Dcreate_vers
#endif
#define H5Dcreate_vers 2

#include <hdf5.h>

#include "bbox.hh"
#include "defs.hh"

#include "tensor_component.hh"



namespace CarpetIOF5 {
  
  namespace F5 {
    
    class data_region_t {
      
      tensor_component_t & m_tensor_component;
      
      bbox<int, dim> const m_region;
      
      hid_t m_properties;
      hid_t m_dataset;
      hid_t m_dataspace;
      
      data_region_t ();
      data_region_t (data_region_t const &);
      data_region_t operator= (data_region_t const &);
      
    public:
      
      data_region_t (tensor_component_t & tensor_component,
                     bbox<int, dim> const & region);
      
      virtual
      ~ data_region_t ();
      
      static string
      name_from_region (bbox<int, dim> const & region);
      
      tensor_component_t &
      get_tensor_component ()
        const;
      
      void
      write (void const * data,
             int cactus_datatype)
        const;
      
      virtual bool
      invariant ()
        const;
      
    };
    
  } // namespace F5

} // namespace CarpetIOF5

#endif  // #ifndef DATA_REGION_HH
