#include <algorithm>
#include <cassert>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include <hdf5.h>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "defs.hh"

#include "file.hh"
#include "utils.hh"

namespace CarpetIOF5 {

namespace F5 {

using namespace std;

file_t::file_t(cGH const *const cctkGH, string const path,
               string const basename, string const extension,
               bool const do_truncate, bool const is_metafile,
               bool const is_datafile)
    : m_cctkGH(cctkGH), m_path(path), m_basename(basename),
      m_extension(extension), m_is_metafile(is_metafile),
      m_is_datafile(is_datafile) {
  DECLARE_CCTK_PARAMETERS;

  assert(cctkGH);

  // Create file names, creating subdirectories to avoid placing
  // too many files into the same directory
  check(CCTK_CreateDirectory(mode, path.c_str()) >= 0);
  if (is_metafile) {
    m_filename = basename + extension;
  } else {
    int const myproc = CCTK_MyProc(cctkGH);
    m_filename = create_filename(myproc, true);
  }
  string const filepath = m_path + "/" + m_filename;

  htri_t is_hdf5;
  H5E_BEGIN_TRY { is_hdf5 = H5Fis_hdf5(filepath.c_str()); }
  H5E_END_TRY;
  bool const file_exists = is_hdf5 > 0;

  if (do_truncate or not file_exists) {
    if (veryverbose) {
      CCTK_VInfo(CCTK_THORNSTRING, "H5Fcreate (name=\"%s\")", filepath.c_str());
    }
    m_hdf5_file =
        H5Fcreate(filepath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  } else {
    if (veryverbose) {
      CCTK_VInfo(CCTK_THORNSTRING, "H5Fopen (name=\"%s\")", filepath.c_str());
    }
    m_hdf5_file = H5Fopen(filepath.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  }
  assert(m_hdf5_file >= 0);

  m_hdf5_fiber_contents = -1;
  m_hdf5_fiber_global_charts = -1;
  m_hdf5_fiber_parameter_space = -1;
  if (m_is_metafile) {
    create_or_check_version(m_hdf5_file);
    m_hdf5_fiber_contents =
        open_or_create_group(m_hdf5_file, "TableOfContents");
    assert(m_hdf5_fiber_contents >= 0);
    m_hdf5_fiber_global_charts = open_or_create_group(m_hdf5_file, "Charts");
    assert(m_hdf5_fiber_global_charts >= 0);
    m_hdf5_fiber_parameter_space =
        open_or_create_group(m_hdf5_fiber_contents, "Parameters");
    assert(m_hdf5_fiber_parameter_space >= 0);

    hid_t const hdf5_time_group =
        open_or_create_group(m_hdf5_fiber_parameter_space, "Time");
    assert(hdf5_time_group >= 0);
    check(not H5Gclose(hdf5_time_group));
  }

  assert(invariant());
}

file_t::~file_t() {
  if (m_is_metafile) {
    check(not H5Gclose(m_hdf5_fiber_parameter_space));
    check(not H5Gclose(m_hdf5_fiber_contents));
    check(not H5Gclose(m_hdf5_fiber_global_charts));
  }
  check(not H5Fclose(m_hdf5_file));
}

string file_t::create_filename(int const proc,
                               bool const create_directories) const {
  DECLARE_CCTK_PARAMETERS;

  int const nprocs = CCTK_nProcs(m_cctkGH);
  string extrapath;
  if (create_subdirs) {
    if (nprocs >= 10000) {
      ostringstream buf;
      buf << extrapath << m_basename << ".p"
          << setw(max(0, processor_digits - 4)) << setfill('0') << proc / 10000
          << "nnnn/";
      extrapath = buf.str();
      if (create_directories) {
        string const filepath = m_path + "/" + extrapath;
        check(CCTK_CreateDirectory(mode, filepath.c_str()) >= 0);
      }
    }
    if (nprocs >= 100) {
      ostringstream buf;
      buf << extrapath << m_basename << ".p"
          << setw(max(0, processor_digits - 2)) << setfill('0') << proc / 100
          << "nn/";
      extrapath = buf.str();
      if (create_directories) {
        string const filepath = m_path + "/" + extrapath;
        check(CCTK_CreateDirectory(mode, filepath.c_str()) >= 0);
      }
    }
  }
  if (one_dir_per_file and nprocs >= 1) {
    ostringstream buf;
    buf << extrapath << m_basename << ".p" << setw(processor_digits)
        << setfill('0') << proc << "/";
    extrapath = buf.str();
    if (create_directories) {
      string const filepath = m_path + "/" + extrapath;
      check(CCTK_CreateDirectory(mode, filepath.c_str()) >= 0);
    }
  }
  ostringstream buf;
  buf << extrapath << m_basename << ".p" << setw(processor_digits)
      << setfill('0') << proc << m_extension;
  return buf.str();
}

void file_t::create_or_check_version(hid_t const hdf5_file) {
  hid_t const hdf5_group = open_or_create_group(hdf5_file, "version");
  assert(hdf5_group >= 0);
  string const version("http://www.zib.de/visual/F5-0.1.2/");
  write_or_check_attribute(hdf5_group, "version", version.c_str());
  check(not H5Gclose(hdf5_group));
}

cGH const *file_t::get_cctkGH() const { return m_cctkGH; }

hid_t file_t::get_hdf5_file() const { return m_hdf5_file; }

hid_t file_t::get_hdf5_fiber_parameter_space() const {
  return m_hdf5_fiber_parameter_space;
}

bool file_t::get_is_metafile() const { return m_is_metafile; }

bool file_t::get_is_datafile() const { return m_is_datafile; }

void file_t::get_link_destination(int const proc, string &filename,
                                  string &objectname) const {
  filename = create_filename(proc);
  objectname = "";
}

bool file_t::invariant() const { return m_cctkGH != 0 and m_hdf5_file >= 0; }

} // namespace F5

} // namespace CarpetIOF5
