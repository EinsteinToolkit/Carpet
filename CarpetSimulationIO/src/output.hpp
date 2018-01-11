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

enum class data_handling { write, attach };

class output_file_t {
  const cGH *const cctkGH;
  const io_dir_t io_dir;
  const string projectname;
  const file_type output_type;
  const int myioproc, ioproc_every;
  shared_ptr<Project> project;
  vector<function<void()> > tasks;

public:
  output_file_t() = delete;
  output_file_t(const output_file_t &) = delete;
  output_file_t(output_file_t &&) = delete;
  output_file_t &operator=(const output_file_t &) = delete;
  output_file_t &operator=(output_file_t &&) = delete;

  output_file_t(const cGH *cctkGH, io_dir_t io_dir, const string &projectname,
                file_type output_type, int myioproc, int ioproc_every);
  ~output_file_t();

  // Insert variables into project
  void insert_vars(const vector<int> &varindices, int reflevel, int timelevel,
                   data_handling handle_data);

  // Write project to file
  void write();
};
} // namespace CarpetSimulationIO

#endif // #ifndef OUTPUT_HPP
