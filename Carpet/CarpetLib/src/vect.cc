// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/vect.cc,v 1.13 2004/01/25 14:57:30 schnetter Exp $

#include <assert.h>

#include <iostream>

#include "defs.hh"

#include "vect.hh"

using namespace std;



// Input
template<class T,int D>
void vect<T,D>::input (istream& is) {
  skipws (is);
  consume (is, '[');
  for (int d=0; d<D; ++d) {
    is >> (*this)[d];
    if (d<D-1) {
      skipws (is);
      consume (is, ',');
    }
  }
  skipws (is);
  consume (is, ']');
}



// Output
template<class T,int D>
void vect<T,D>::output (ostream& os) const {
  os << "[";
  for (int d=0; d<D; ++d) {
    os << (*this)[d];
    if (d<D-1) os << ",";
  }
  os << "]";
}



// Note: We need all dimensions all the time.
template class vect<int,1>;
template class vect<int,2>;
template class vect<int,3>;

template void vect<double,3>::input (istream& is);
template void vect<vect<bool,2>,3>::input (istream& is);

template void vect<bool,2>::output (ostream& os) const;
template void vect<bool,3>::output (ostream& os) const;
template void vect<double,3>::output (ostream& os) const;
template void vect<vect<bool,2>,3>::output (ostream& os) const;
template void vect<vect<int,2>,3>::output (ostream& os) const;
