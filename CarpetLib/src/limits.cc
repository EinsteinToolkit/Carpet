#include <cctk.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <iostream>
#ifdef HAVE_SYS_RESOURCE_H
#include <sys/resource.h>
#endif

#include "defs.hh"

#include "limits.hh"

namespace CarpetLib {
using namespace std;

#ifdef HAVE_SYS_RESOURCE_H

static void set_limit(int resource, char const *name, CCTK_INT value);

static ostream &operator<<(ostream &s, struct rlimit const &limit);

static void output(ostream &s, rlim_t const &value);

void set_system_limits() {
  DECLARE_CCTK_PARAMETERS;
  set_limit(RLIMIT_CORE, "core file size", max_core_size_MB);
  set_limit(RLIMIT_AS, "addres space size", max_memory_size_MB);
  set_limit(RLIMIT_DATA, "data segment size", max_memory_size_MB);
#ifdef RLIMIT_RSS
  set_limit(RLIMIT_RSS, "resident set size", max_memory_size_MB);
#endif
}

void set_limit(int const resource, char const *const name,
               CCTK_INT const value) {
  int ierr;
  struct rlimit limit;
  ierr = getrlimit(resource, &limit);
  if (ierr < 0)
    CCTK_ERROR("getrlimit failed");

  if (value == -2) {
    // Only show limit
    cout << "Current " << name << " limit: " << limit << endl;
    return;
  }

  cout << "Old " << name << " limit: " << limit << endl;

  if (value == -1) {
    limit.rlim_cur = limit.rlim_max;
  } else {
    limit.rlim_cur = min((rlim_t)value * 1024 * 1024, limit.rlim_max);
  }

  ierr = setrlimit(resource, &limit);
  if (ierr < 0)
    CCTK_ERROR("setrlimit failed");
  ierr = getrlimit(resource, &limit);
  if (ierr < 0)
    CCTK_ERROR("getrlimit failed");

  cout << "New " << name << " limit: " << limit << endl;
}

static ostream &operator<<(ostream &s, struct rlimit const &limit) {
  s << "hard=";
  output(s, limit.rlim_max);
  s << ", soft=";
  output(s, limit.rlim_cur);
  return s;
}

static void output(ostream &s, rlim_t const &value) {
  if (value == RLIM_INFINITY) {
    s << "[unlimited]";
  } else {
    s << (value / CCTK_REAL(1024 * 1024)) << " MB";
  }
}

#else

void set_system_limits() { return; }

void set_limit(int const /*resource*/, char const *const /*name*/,
               CCTK_INT const /*value*/) {
  return;
}

#endif

} // namespace CarpetLib
