#ifndef OUTPUT_HPP
#define OUTPUT_HPP

#include "util.hpp"

#include <cctk.h>

#include <SimulationIO/SimulationIO.hpp>

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace CarpetSimulationIO {
using namespace SimulationIO;
using namespace std;

class output_file_t {
  const cGH *cctkGH;
  io_dir_t io_dir;
  string projectname;
  int ioproc, nioprocs;
  shared_ptr<Project> project;
  vector<function<void()> > tasks;

public:
  output_file_t() = delete;
  output_file_t(const output_file_t &) = delete;
  output_file_t(output_file_t &&) = delete;
  output_file_t &operator=(const output_file_t &) = delete;
  output_file_t &operator=(output_file_t &&) = delete;

  output_file_t(const cGH *cctkGH, io_dir_t io_dir, const string &projectname,
                int ioproc, int nioprocs);
  ~output_file_t();

  // Insert variables into project
  void insert_vars(const vector<int> &varindices, int reflevel, int timelevel,
                   file_type outfile);

  // Write project to file
  void write();
};
} // namespace CarpetSimulationIO

#endif // #ifndef OUTPUT_HPP
