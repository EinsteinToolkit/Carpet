#ifndef UTIL_HH
#define UTIL_HH

#include <iostream>
#include <vector>

namespace Requirements {

using namespace std;

  // taken from defs.cc and defs.hh
  // Vector output
  template<class T>
  inline ostream& output (ostream& os, const vector<T>& v) {
    os << "[";
    // Do not number the elements, as this would lead to a format that
    // cannot be read back in.
  //   int cnt=0;
    for (typename vector<T>::const_iterator ti=v.begin(); ti!=v.end(); ++ti) {
      if (ti!=v.begin()) os << ",";
  //     os << cnt++ << ":";
      os << *ti;
    }
    os << "]";
    return os;
  }

  template<class T>
  inline ostream& operator<< (ostream& os, const vector<T>& v) {
    return Requirements::output(os,v);
  }

};

#endif
