#ifndef LOCATION_HH
#define LOCATION_HH

#include <iostream>

namespace Requirements {

  using namespace std;

  // Struct defining a location of a grid point
  struct location_t {
    int it, vi, tl, rl, m;
    char const* info;
    location_t():
      it(-1), vi(-1), tl(-1), rl(-1), m(-1), info("")
    {}
    location_t(int _it, int _vi, int _tl, int _rl, int _m, char const* _info):
      it(_it), vi(_vi), tl(_tl), rl(_rl), m(_m), info(_info)
    {}
    // Output helper
    void output (ostream& os) const;
  };

  inline ostream& operator<< (ostream& os, const location_t& a) {
    a.output(os);
    return os;
  }
}

#endif
