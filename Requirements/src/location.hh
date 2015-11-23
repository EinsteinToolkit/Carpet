#ifndef LOCATION_HH
#define LOCATION_HH

#include <cctk.h>
#include <cctk_Schedule.h>

#include <iostream>
#include <string>

namespace Requirements {

using namespace std;

// Struct defining a location of a grid point
struct location_t {
  string info;
  cFunctionData const *fd;
  int it;
  int vi;
  int rl, m, tl;
  location_t() : info(""), fd(0), it(-1), vi(-1), rl(-1), m(-1), tl(-1) {}
  location_t(string info_)
      : info(info_), fd(0), it(-1), vi(-1), rl(-1), m(-1), tl(-1) {}
  location_t(string info_, cFunctionData const *fd_)
      : info(info_), fd(fd_), it(-1), vi(-1), rl(-1), m(-1), tl(-1) {}
  // Output helper
  void output(ostream &os) const;
};

inline ostream &operator<<(ostream &os, const location_t &loc) {
  loc.output(os);
  return os;
}
}

#endif
