#ifndef INPUT_HPP
#define INPUT_HPP

#include "util.hpp"

#include <cctk.h>

#include <SimulationIO/SimulationIO.hpp>

#include <memory>
#include <string>
#include <vector>

namespace CarpetSimulationIO {
using namespace SimulationIO;
using namespace std;

class input_file_t {
  const cGH *const cctkGH;
  const int iteration;
  shared_ptr<Project> project;

public:
  input_file_t() = delete;
  input_file_t(const input_file_t &) = delete;
  input_file_t(input_file_t &&) = delete;
  input_file_t &operator=(const input_file_t &) = delete;
  input_file_t &operator=(input_file_t &&) = delete;

  input_file_t(const cGH *cctkGH, io_dir_t io_dir, const string &projectname,
               int iteration, int ioproc, int nioprocs);

  // Read parameters from project
  void read_params() const;

  // Read variables from project
  void read_vars(const vector<int> &varindices, int reflevel,
                 int timelevel) const;
};
} // namespace CarpetSimulationIO

#endif // #ifndef INPUT_HPP
