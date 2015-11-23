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

  // File mode for creating directories
  static int const mode = 0755;

  cGH const *const m_cctkGH;

  string const m_path;
  string const m_basename;
  string const m_extension;
  string m_filename;

  bool const m_is_metafile;
  bool const m_is_datafile;

  hid_t m_hdf5_file;

  hid_t m_hdf5_fiber_contents;
  hid_t m_hdf5_fiber_global_charts;
  hid_t m_hdf5_fiber_parameter_space;

  file_t();
  file_t(file_t const &);
  file_t operator=(file_t const &);

  string create_filename(int proc, bool create_directories = false) const;

  static void create_or_check_version(hid_t const hdf5_file);

public:
  file_t(cGH const *cctkGH, string path, string basename, string extension,
         bool do_truncate, bool is_metafile, bool is_datafile);

  virtual ~file_t();

  cGH const *get_cctkGH() const;

  hid_t get_hdf5_file() const;

  hid_t get_hdf5_fiber_parameter_space() const;

  bool get_is_metafile() const;

  bool get_is_datafile() const;

  void get_link_destination(int proc, string &filename,
                            string &objectname) const;

  virtual bool invariant() const;
};

} // namespace F5

} // namespace CarpetIOF5

#endif // #ifndef FILE_HH
