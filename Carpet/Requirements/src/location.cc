#include <iostream>

#include <location.hh>

namespace Requirements {

  void location_t::output(ostream& os) const
  {
    os << "vi:" << vi << ", "
       << "["
       << "it:" << it << ","
       << "rl:" << rl << ","
       << "m:"  << m  << ","
       << "tl:" << tl
       << "]";
  }
}
