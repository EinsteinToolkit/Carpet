#ifndef TIMESTEP_HH
#define TIMESTEP_HH

#include <string>

#include <hdf5.h>

#include "cctk.h"

#include "file.hh"

namespace CarpetIOF5 {

using std::string;

namespace F5 {

class timestep_t {

  file_t &m_file;

  CCTK_REAL const m_time;

  string m_name;

  hid_t m_hdf5_timestep;

  timestep_t();
  timestep_t(timestep_t const &);
  timestep_t operator=(timestep_t const &);

public:
  timestep_t(file_t &file, CCTK_REAL time, char const *name = 0);

  virtual ~timestep_t();

  file_t &get_file() const;

  CCTK_REAL
  get_time() const;

  hid_t get_hdf5_timestep() const;

  void get_link_destination(int proc, string &filename,
                            string &objectname) const;

  virtual bool invariant() const;
};

} // namespace F5

} // namespace CarpetIOF5

#endif // #ifndef TIMESTEP_HH
