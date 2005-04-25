#include <cassert>

#include <hdf5.h>

#include "cctk.h"

#include "file.hh"



namespace CarpetIOF5 {
  
  namespace F5 {
    
    file_t::
    file_t (cGH const * const cctkGH,
            char const * const filename)
      : m_cctkGH (cctkGH),
        m_filename (filename)
    {
      assert (cctkGH);
      assert (filename);
      
      extending extension (cctkGH);
      
      bool const should_truncate = IO_TruncateOutputFiles (cctkGH);
      
      bool const did_truncate
        = (extension.get_did_truncate().find (filename)
           != extension.get_did_truncate().end());
      
      htri_t is_hdf5 = H5Fis_hdf5 (filename);
      assert (is_hdf5 >= 0);
      bool const file_exists = is_hdf5 > 0;
      
      if ((should_truncate and ! did_truncate) or ! file_exists)
      {
        m_hdf5_file
          = H5Fcreate (filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        extension.get_did_truncate().insert (filename);
      }
      else
      {
        m_hdf5_file = H5Fopen (filename, H5F_ACC_RDWR, H5P_DEFAULT);
      }
      
      assert (invariant());
    }
    
    
    
    file_t::
    ~ file_t()
    {
      herr_t const herr = H5Fclose (m_hdf5_file);
      assert (! herr);
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
              and m_filename != 0
              and m_hdf5_file >= 0);
    }
    
  } // namespace F5

} // namespace CarpetIOF5
