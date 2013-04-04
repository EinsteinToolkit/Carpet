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
    location_t(int it_, int vi_, int tl_, int rl_, int m_, char const* info_):
      it(it_), vi(vi_), tl(tl_), rl(rl_), m(m_), info(info_)
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
