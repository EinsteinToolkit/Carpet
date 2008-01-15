#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <hdf5.h>

#include "cctk.h"

#include "file.hh"



namespace CarpetIOF5 {
  
  namespace F5 {
    
    using namespace std;
    
    file_t::
    file_t (cGH const * const cctkGH,
            string const basename,
            string const extension,
            bool const want_metafile,
            bool const do_truncate)
      : m_cctkGH (cctkGH),
        m_have_metafile (want_metafile),
        m_basename (basename),
        m_extension (extension)
    {
      assert (cctkGH);
      
      int const proc = CCTK_MyProc (cctkGH);
      m_metafilename = make_metafilename ();
      m_filename = make_filename (proc);
      
      char const * const metafilenameptr = m_metafilename.c_str();
      char const * const filenameptr = m_filename.c_str();
      
      htri_t is_hdf5;
      H5E_BEGIN_TRY {
        is_hdf5 = H5Fis_hdf5 (filenameptr);
      } H5E_END_TRY;
      bool const file_exists = is_hdf5 > 0;
      
      // Ignore the metafile when determining whether the file already
      // exists
      if (do_truncate or not file_exists)
      {
        m_hdf5_metafile = -1;
        if (m_have_metafile)
        {
          m_hdf5_metafile
            = H5Fcreate (metafilenameptr,
                         H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        }
        m_hdf5_file
          = H5Fcreate (filenameptr, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      }
      else
      {
        m_hdf5_metafile = -1;
        if (m_have_metafile)
        {
          m_hdf5_metafile
            = H5Fopen (metafilenameptr, H5F_ACC_RDWR, H5P_DEFAULT);
        }
        m_hdf5_file = H5Fopen (filenameptr, H5F_ACC_RDWR, H5P_DEFAULT);
      }
      
      m_filenames.resize (CCTK_nProcs (cctkGH));
      m_filenames.at(proc) = m_filename;
      
      assert (invariant());
    }
    
    
    
    file_t::
    ~ file_t()
    {
      if (m_have_metafile)
      {
        herr_t const herr = H5Fclose (m_hdf5_metafile);
        assert (not herr);
      }
      herr_t const herr = H5Fclose (m_hdf5_file);
      assert (not herr);
    }
    
    
    
    int file_t::
    base_10_digits (int number)
    {
      number = abs (number);
      int digits = 1;
      while (number >= 10)
      {
        number /= 10;
        ++ digits;
      }
      return digits;
    }
    
    
    
    string file_t::
    make_metafilename ()
      const
    {
      return m_basename + m_extension;
    }
    
    
    
    string file_t::
    make_filename (int const proc)
      const
    {
      int const digits = base_10_digits (CCTK_nProcs (m_cctkGH) - 1);
      ostringstream filenamebuf;
      filenamebuf << m_basename
                  << "."
                  << setw (digits) << setfill ('0') << proc
                  << m_extension;
      return filenamebuf.str();
    }
    
    
    
    cGH const * file_t::
    get_cctkGH ()
      const
    {
      return m_cctkGH;
    }
    
    
    
    bool file_t::
    get_have_metafile ()
      const
    {
      return m_have_metafile;
    }
    
    
    
    string file_t::
    get_filename (int const proc)
      const
    {
      assert (proc >= 0 and proc < CCTK_nProcs (m_cctkGH));
      if (m_filenames.at(proc).empty())
      {
        m_filenames.at(proc) = make_filename (proc);
      }
      return m_filenames.at(proc);
    }
    
    
    
    hid_t file_t::
    get_hdf5_metafile()
      const
    {
      return m_hdf5_metafile;
    }
    
    
    
    hid_t file_t::
    get_hdf5_file()
      const
    {
      return m_hdf5_file;
    }
    
    
    
    void file_t::
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
        l_filename = m_filename;
        l_objectname = string ("");
      }
      filename = l_filename;
      objectname = l_objectname;
    }
    
    
    
    bool file_t::
    invariant()
      const
    {
      return (m_cctkGH != 0
              and (not m_have_metafile or m_hdf5_metafile >= 0)
              and m_hdf5_file >= 0);
    }
    
  } // namespace F5

} // namespace CarpetIOF5
