#include "location.hh"

#include <cctk.h>

#include <cstdlib>
#include <iostream>

namespace Requirements {

void location_t::output(ostream &os) const {
  os << "LOC: " << info << " ";
  if (fd) {
    os << "func " << fd->thorn << "::" << fd->routine << " "
       << "in " << fd->where << " ";
  }
  char *const fullname = CCTK_FullName(vi);
  os << "it " << it << ", var " << fullname << " "
     << "["
     << "rl:" << rl << ","
     << "m:" << m << ","
     << "tl:" << tl << "]";
  free(fullname);
}
}
