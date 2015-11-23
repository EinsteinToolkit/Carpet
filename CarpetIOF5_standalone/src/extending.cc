#include <cassert>
#include <cstring>

#include "cctk.h"

#include "extending.hh"

namespace CarpetIOF5 {

char const *const extending_t::extension_name = "CarpetIOF5";

int extending_t::create(void *(*const Setup)(tFleshConfig *config,
                                             int convlevel, cGH *cctkGH)) {
  int const handle = CCTK_RegisterGHExtension(extension_name);
  assert(handle >= 0);
  int const iflag = CCTK_RegisterGHExtensionSetupGH(handle, Setup);
  assert(iflag);
  return 0; // no error
}

void *extending_t::setup(
    cGH *const cctkGH, int (*const OutputGH)(cGH const *cctkGH),
    int (*const TimeToOutput)(cGH const *cctkGH, int variable),
    int (*const TriggerOutput)(cGH const *cctkGH, int variable),
    int (*const OutputVarAs)(cGH const *cctkGH, char const *varname,
                             char const *alias)) {
  assert(cctkGH != 0);

  int ierr;
  int const io_method = CCTK_RegisterIOMethod(extension_name);
  ierr = CCTK_RegisterIOMethodOutputGH(io_method, OutputGH);
  assert(not ierr);
  ierr = CCTK_RegisterIOMethodTimeToOutput(io_method, TimeToOutput);
  assert(not ierr);
  ierr = CCTK_RegisterIOMethodTriggerOutput(io_method, TriggerOutput);
  assert(not ierr);
  ierr = CCTK_RegisterIOMethodOutputVarAs(io_method, OutputVarAs);
  assert(not ierr);

  return new extension_t;
}

extending_t::extending_t(cGH const *cctkGH) {
  assert(cctkGH);
  void *const ext = CCTK_GHExtension(cctkGH, extension_name);
  assert(ext != 0);
  m_extension = static_cast<extension_t *>(ext);
}

bool extending_t::get_did_truncate(string const name) const {
  return (m_extension->did_truncate.find(name) !=
          m_extension->did_truncate.end());
}

void extending_t::set_did_truncate(string const name) {
  m_extension->did_truncate.insert(name);
}

int extending_t::get_last_output_iteration(int const ml, int const rl,
                                           int const vi) const {
  resize_last_output_iteration(ml, rl, vi);
  return m_extension->last_output_iteration.at(ml).at(rl).at(vi);
}

void extending_t::set_last_output_iteration(int const ml, int const rl,
                                            int const vi, int const iteration) {
  resize_last_output_iteration(ml, rl, vi);
  m_extension->last_output_iteration.at(ml).at(rl).at(vi) = iteration;
}

void extending_t::resize_last_output_iteration(int ml, int rl, int vi) const {
  assert(ml >= 0);
  if (size_t(ml) >= m_extension->last_output_iteration.size()) {
    m_extension->last_output_iteration.resize(ml + 1);
  }
  assert(rl >= 0);
  if (size_t(rl) >= m_extension->last_output_iteration.at(ml).size()) {
    m_extension->last_output_iteration.at(ml).resize(rl + 1);
  }
  if (size_t(vi) >= m_extension->last_output_iteration.at(ml).at(rl).size()) {
    m_extension->last_output_iteration.at(ml).at(rl).resize(vi + 1, -1);
  }
}

} // namespace CarpetIOF5
