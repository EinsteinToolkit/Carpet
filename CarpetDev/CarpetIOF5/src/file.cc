#include <cassert>
#include <string>

#include <hdf5.h>

#include "cctk.h"

#include "file.hh"



namespace CarpetIOF5 {
  
  namespace F5 {
    
    file_t::
    file_t (cGH const * const cctkGH,
            string const filename,
            bool const do_truncate)
      : m_cctkGH (cctkGH),
        m_filename (filename)
    {
      assert (cctkGH);
      
      char const * const filenameptr = filename.c_str();
      
      htri_t is_hdf5;
      H5E_BEGIN_TRY {
        is_hdf5 = H5Fis_hdf5 (filenameptr);
      } H5E_END_TRY;
      bool const file_exists = is_hdf5 > 0;
      
      if (do_truncate or not file_exists)
      {
        m_hdf5_file
          = H5Fcreate (filenameptr, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      }
      else
      {
        m_hdf5_file = H5Fopen (filenameptr, H5F_ACC_RDWR, H5P_DEFAULT);
      }
      
      assert (invariant());
    }
    
    
    
    file_t::
    ~ file_t()
    {
      herr_t const herr = H5Fclose (m_hdf5_file);
      assert (not herr);
    }
    
    
    
    cGH const * file_t::
    get_cctkGH ()
      const
    {
      return m_cctkGH;
    }
    
    
    
    hid_t file_t::
    get_hdf5_file()
      const
    {
      return m_hdf5_file;
    }
    
    
    
    bool file_t::
    invariant()
      const
    {
      return (m_cctkGH != 0
              and m_hdf5_file >= 0);
    }
    
  } // namespace F5

} // namespace CarpetIOF5
