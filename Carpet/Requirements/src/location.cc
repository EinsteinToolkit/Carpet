#include <location.hh>
#include <iostream>

namespace Requirements {

  void location_t::output(ostream& os) const
  {
    os << "vi:"  << vi << ",it:" << it << ", ";
    os << "[rl:" << rl << ",";
    os <<  "tl:" << tl << ",";
    os <<  "m:"  << m << "]";
  }
}
