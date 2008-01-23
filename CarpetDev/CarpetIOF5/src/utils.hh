#ifndef UTILS_HH
#define UTILS_HH

// force HDF5 1.8.x installations to use the new API
#ifdef H5Gopen_vers
#undef H5Gopen_vers
#endif
#ifdef H5Gcreate_vers
#undef H5Gcreate_vers
#endif
#ifdef H5Acreate_vers
#undef H5Acreate_vers
#endif
#ifdef H5Tarray_create_vers
#undef H5Tarray_create_vers
#endif
#define H5Gopen_vers 2
#define H5Gcreate_vers 2
#define H5Acreate_vers 2
#define H5Tarray_create_vers 2

#include <hdf5.h>

#include "cctk.h"

#include "vect.hh"



namespace CarpetIOF5 {
  
  namespace F5 {
    
    hid_t
    hdf5_datatype_from_dummy (signed char const & dummy);
    hid_t
    hdf5_datatype_from_dummy (short const & dummy);
    hid_t
    hdf5_datatype_from_dummy (int const & dummy);
    hid_t
    hdf5_datatype_from_dummy (long const & dummy);
    hid_t
    hdf5_datatype_from_dummy (long long const & dummy);
    hid_t
    hdf5_datatype_from_dummy (float const & dummy);
    hid_t
    hdf5_datatype_from_dummy (double const & dummy);
    hid_t
    hdf5_datatype_from_dummy (long double const & dummy);
#ifdef HAVE_CCTK_COMPLEX8
    hid_t
    hdf5_datatype_from_dummy (CCTK_COMPLEX8 const & dummy);
#endif
#ifdef HAVE_CCTK_COMPLEX16
    hid_t
    hdf5_datatype_from_dummy (CCTK_COMPLEX16 const & dummy);
#endif
#ifdef HAVE_CCTK_COMPLEX32
    hid_t
    hdf5_datatype_from_dummy (CCTK_COMPLEX32 const & dummy);
#endif
    
    template<typename T, typename R>
    hid_t
    hdf5_complex_datatype_from_dummy (T const & dummy, R const & real);
    
    hid_t
    hdf5_datatype_from_cactus_datatype (int cactus_datatype);
    
    
    
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
