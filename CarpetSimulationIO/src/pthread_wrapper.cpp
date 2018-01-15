#include "pthread_wrapper.hpp"

#include <cctk_Parameters.h>

#include <cassert>
#include <iostream>

namespace CarpetSimulationIO {
using namespace std;

#if HAVE_CAPABILITY_PTHREADS

void *pthread_wrapper_t::run(void *arg) {
  return static_cast<pthread_wrapper_t *>(arg)->run1();
}

void *pthread_wrapper_t::run1() {
  DECLARE_CCTK_PARAMETERS;
  if (verbose)
    CCTK_VINFO("Starting thread \"%s\"", name.c_str());
  fun();
  if (verbose)
    CCTK_VINFO("Finished thread \"%s\"", name.c_str());
  return nullptr;
}

pthread_wrapper_t::pthread_wrapper_t(const string &name, const fun_t &fun)
    : name(name), fun(fun) {
  DECLARE_CCTK_PARAMETERS;
  if (verbose)
    CCTK_VINFO("Creating thread \"%s\"", name.c_str());
  int ierr = pthread_create(&thread, NULL, &run, this);
  assert(not ierr);
}

pthread_wrapper_t::~pthread_wrapper_t() {
  DECLARE_CCTK_PARAMETERS;
  if (verbose)
    CCTK_VINFO("Waiting for thread \"%s\"", name.c_str());
  void *retval;
  int ierr = pthread_join(thread, &retval);
  assert(not ierr);
}

#endif

} // namespace CarpetSimulationIO
