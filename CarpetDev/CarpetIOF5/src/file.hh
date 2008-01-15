#ifndef FILE_HH
#define FILE_HH

#include <string>
#include <vector>

#include <hdf5.h>

#include "cctk.h"

#include "extending.hh"



namespace CarpetIOF5 {
  
  using std::string;
  using std::vector;
  
  
  
  namespace F5 {
    
    class file_t {
      
      cGH const * const m_cctkGH;
      
      bool const m_have_metafile;
      
      string const m_basename;
      string const m_extension;
      
      string m_metafilename;
      string m_filename;
      
      mutable vector <string> m_filenames;
      
      hid_t m_hdf5_metafile;
      hid_t m_hdf5_file;
      
      file_t ();
      file_t (file_t const &);
      file_t operator= (file_t const &);
      
      static int
      base_10_digits (int number);
      
      string
      make_metafilename ()
        const;
      
      string
      make_filename (int proc)
        const;
      
    public:
      
      file_t (cGH const * cctkGH,
              string filename,
              string extension,
              bool want_metafile,
              bool do_truncate);
      
      virtual
      ~ file_t ();
      
      cGH const *
      get_cctkGH ()
        const;
      
      bool
      get_have_metafile ()
        const;
      
      string
      get_filename (int proc)
        const;
      
      hid_t
      get_hdf5_metafile ()
        const;
      
      hid_t
      get_hdf5_file ()
        const;
      
      void
      get_link_destination (string & filename,
                            string & objectname)
        const;
      
      virtual bool
      invariant ()
        const;
    };
    
  } // namespace F5

} // namespace CarpetIOF5

#endif  // #ifndef FILE_HH
