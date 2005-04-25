#include <cassert>
#include <complex>
#include <vector>

#include <hdf5.h>

#include "cctk.h"

#include "defs.hh"

#include "utils.hh"



namespace CarpetIOF5 {
  
  namespace F5 {
    
    using std::complex;
    using std::vector;
    
    
    
    hid_t
    hdf5_datatype (signed char const & dummy)
    {
      return H5T_NATIVE_SCHAR;
    }
    
    hid_t
    hdf5_datatype (short const & dummy)
    {
      return H5T_NATIVE_SHORT;
    }
    
    hid_t
    hdf5_datatype (int const & dummy)
    {
      return H5T_NATIVE_INT;
    }
    
    hid_t
    hdf5_datatype (long const & dummy)
    {
      return H5T_NATIVE_LONG;
    }
    
    hid_t
    hdf5_datatype (long long const & dummy)
    {
      return H5T_NATIVE_LLONG;
    }
    
    hid_t
    hdf5_datatype (float const & dummy)
    {
      return H5T_NATIVE_FLOAT;
    }
    
    hid_t
    hdf5_datatype (double const & dummy)
    {
      return H5T_NATIVE_DOUBLE;
    }
    
    hid_t
    hdf5_datatype (long double const & dummy)
    {
      return H5T_NATIVE_LDOUBLE;
    }
    
    hid_t
    hdf5_datatype (CCTK_COMPLEX8 const & dummy)
    {
      CCTK_REAL4 real;
      return hdf5_complex_datatype (dummy, real);
    }
    
    hid_t
    hdf5_datatype (CCTK_COMPLEX16 const & dummy)
    {
      CCTK_REAL8 real;
      return hdf5_complex_datatype (dummy, real);
    }
    
    hid_t
    hdf5_datatype (CCTK_COMPLEX32 const & dummy)
    {
      CCTK_REAL16 real;
      return hdf5_complex_datatype (dummy, real);
    }
    
    template<typename T, typename R>
    hid_t
    hdf5_complex_datatype (T const & dummy, R const & real)
    {
      static bool initialised = false;
      static hid_t hdf_complex;
      
      if (! initialised)
      {
        initialised = true;
        
        hsize_t const dim = 2;
        int const perm = 0;
        
        hdf_complex = H5Tarray_create (hdf5_datatype (real), 1, & dim, & perm);
        assert (hdf_complex >= 0);
      }
      
      return hdf_complex;
    }
    
    
    
    hid_t
    hdf5_datatype_from_cactus (int const cactus_type)
    {
      switch (cactus_type) {
      case CCTK_VARIABLE_VOID     : return H5I_INVALID_HID;
      case CCTK_VARIABLE_BYTE     : { CCTK_BYTE      const dummy = 0; return hdf5_datatype (dummy); }
      case CCTK_VARIABLE_INT      : { CCTK_INT       const dummy = 0; return hdf5_datatype (dummy); }
      case CCTK_VARIABLE_INT1     : { CCTK_INT1      const dummy = 0; return hdf5_datatype (dummy); }
      case CCTK_VARIABLE_INT2     : { CCTK_INT2      const dummy = 0; return hdf5_datatype (dummy); }
      case CCTK_VARIABLE_INT4     : { CCTK_INT4      const dummy = 0; return hdf5_datatype (dummy); }
      case CCTK_VARIABLE_INT8     : { CCTK_INT8      const dummy = 0; return hdf5_datatype (dummy); }
      case CCTK_VARIABLE_REAL     : { CCTK_REAL      const dummy = 0; return hdf5_datatype (dummy); }
      case CCTK_VARIABLE_REAL4    : { CCTK_REAL4     const dummy = 0; return hdf5_datatype (dummy); }
      case CCTK_VARIABLE_REAL8    : { CCTK_REAL8     const dummy = 0; return hdf5_datatype (dummy); }
      case CCTK_VARIABLE_REAL16   : { CCTK_REAL16    const dummy = 0; return hdf5_datatype (dummy); }
      case CCTK_VARIABLE_COMPLEX  : { CCTK_COMPLEX         dummy    ; return hdf5_datatype (dummy); }
      case CCTK_VARIABLE_COMPLEX8 : { CCTK_COMPLEX8        dummy    ; return hdf5_datatype (dummy); }
      case CCTK_VARIABLE_COMPLEX16: { CCTK_COMPLEX16       dummy    ; return hdf5_datatype (dummy); }
      case CCTK_VARIABLE_COMPLEX32: { CCTK_COMPLEX32       dummy    ; return hdf5_datatype (dummy); }
      case CCTK_VARIABLE_CHAR     : return H5I_INVALID_HID;
      case CCTK_VARIABLE_STRING   : return H5I_INVALID_HID;
      case CCTK_VARIABLE_POINTER  : return H5I_INVALID_HID;
      case CCTK_VARIABLE_FPOINTER : return H5I_INVALID_HID;
      }
      return H5I_INVALID_HID;
    }
    
    
    
    hid_t
    open_or_create_group (hid_t const where,
                          char const * const name)
    {
      assert (where >= 0);
      assert (name != 0);
      
      bool group_exists;
      H5E_BEGIN_TRY {
        H5G_stat_t statbuf;
        herr_t const herr
          = H5Gget_objinfo (where, name, true, & statbuf);
        group_exists = herr == 0;
      } H5E_END_TRY;
      
      hid_t group;
      if (group_exists)
      {
        group = H5Gopen (where, name);
      }
      else
      {
        group = H5Gcreate (where, name, 0);
      }
      
      return group;
    }
    
    
    
    template<typename T>
    void
    write_or_check_attribute (hid_t const where,
                              char const * const name,
                              T const * const values,
                              int const num_values)
    {
      assert (where >= 0);
      assert (name != 0);
      assert (num_values >= 0);
      assert (num_values == 0 or values != 0);
      
      T dummy;
      hid_t const datatype = hdf5_datatype (dummy);
      
      hid_t attribute;
      H5E_BEGIN_TRY {
        attribute = H5Aopen_name (where, name);
      } H5E_END_TRY;
      
      if (attribute < 0)
      {
        // The attribute does not yet exist; create it
        hsize_t const dim = num_values;
        hid_t const dataspace = H5Screate_simple (1, & dim, & dim);
        assert (dataspace >= 0);
        attribute = H5Acreate (where, name, datatype, dataspace, H5P_DEFAULT);
        assert (attribute >= 0);
        herr_t herr;
        herr = H5Awrite (attribute, datatype, values);
        assert (! herr);
        herr = H5Aclose (attribute);
        assert (! herr);
        herr = H5Sclose (dataspace);
        assert (! herr);
      }
      else
      {
        // The attribute already exists; read and check it
        hid_t const dataspace = H5Aget_space (attribute);
        htri_t const is_simple = H5Sis_simple (dataspace);
        assert (is_simple >= 0);
        assert (is_simple > 0);
        int const ndims = H5Sget_simple_extent_ndims (dataspace);
        assert (ndims == 1);
        hsize_t dim;
        herr_t herr;
        herr = H5Sget_simple_extent_dims (dataspace, & dim, 0);
        assert (dim == num_values);
        vector<T> buf (dim);
        herr = H5Aread (attribute, datatype, & buf.front());
        assert (! herr);
        herr = H5Sclose (dataspace);
        assert (! herr);
        herr = H5Aclose (attribute);
        assert (! herr);
        for (int n = 0; n < num_values; ++ n)
        {
          assert (values [n] == buf [n]);
        }
      }
    }
    
    template
    void
    write_or_check_attribute (hid_t const where,
                              char const * const name,
                              int const * const values,
                              int const num_values);
    template
    void
    write_or_check_attribute (hid_t const where,
                              char const * const name,
                              CCTK_REAL const * const values,
                              int const num_values);
    
    
    
    template<typename T>
    void
    write_or_check_attribute (hid_t const where,
                              char const * const name,
                              T const & value)
    {
      assert (where >= 0);
      assert (name != 0);
      
      write_or_check_attribute (where, name, & value, 1);
    }
    
    template
    void
    write_or_check_attribute (hid_t const where,
                              char const * const name,
                              int const & value);
    template
    void
    write_or_check_attribute (hid_t const where,
                              char const * const name,
                              CCTK_REAL const & value);
    
    
    
    template<typename T, int D>
    void
    write_or_check_attribute (hid_t where,
                              char const * name,
                              vect<T,D> const & value)
    {
      assert (where >= 0);
      assert (name != 0);
      
      write_or_check_attribute (where, name, & value [0], D);
    }
    
    template
    void
    write_or_check_attribute (hid_t where,
                              char const * name,
                              vect<int, dim> const & value);
    template
    void
    write_or_check_attribute (hid_t where,
                              char const * name,
                              vect<CCTK_REAL, dim> const & value);

  } // namespace F5

} // namespace CarpetIOF5
