#ifndef utils_HH
#define utils_HH

#include <hdf5.h>

#include "cctk.h"

#include "vect.hh"



namespace CarpetIOF5 {
  
  namespace F5 {
    
    hid_t
    hdf5_datatype (signed char const & dummy);
    hid_t
    hdf5_datatype (short const & dummy);
    hid_t
    hdf5_datatype (int const & dummy);
    hid_t
    hdf5_datatype (long const & dummy);
    hid_t
    hdf5_datatype (long long const & dummy);
    hid_t
    hdf5_datatype (float const & dummy);
    hid_t
    hdf5_datatype (double const & dummy);
    hid_t
    hdf5_datatype (long double const & dummy);
    hid_t
    hdf5_datatype (CCTK_COMPLEX8 const & dummy);
    hid_t
    hdf5_datatype (CCTK_COMPLEX16 const & dummy);
    hid_t
    hdf5_datatype (CCTK_COMPLEX32 const & dummy);
    
    template<typename T, typename R>
    hid_t
    hdf5_complex_datatype (T const & dummy, R const & real);
    
    hid_t
    hdf5_datatype_from_cactus (int cactus_type);
    
    
    
    hid_t
    open_or_create_group (hid_t where,
                          char const * name);
    
    
    
    template<typename T>
    void
    write_or_check_attribute (hid_t where,
                              char const * name,
                              T const * values,
                              int num_values);
    
    template<typename T>
    void
    write_or_check_attribute (hid_t where,
                              char const * name,
                              T const & value);
    
    template<typename T, int D>
    void
    write_or_check_attribute (hid_t where,
                              char const * name,
                              vect<T,D> const & value);
    
  } // namespace F5

} // namespace CarpetIOF5

#endif  // #ifndef UTILS_HH
