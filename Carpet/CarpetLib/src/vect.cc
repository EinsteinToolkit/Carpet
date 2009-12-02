#include <cassert>
#include <iostream>

#include "cctk.h"

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
    assert (is.good());
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



// Specialise some constructors for lower dimensions
// These functions are declared, but must not be used.

template<> vect<int,0>::vect (const int& x, const int& y) { assert(0); }
template<> vect<int,1>::vect (const int& x, const int& y) { assert(0); }
template<> vect<int,3>::vect (const int& x, const int& y) { assert(0); }
template<> vect<int,4>::vect (const int& x, const int& y) { assert(0); }

template<> vect<int,0>::vect (const int& x, const int& y, const int& z) { assert(0); }
template<> vect<int,1>::vect (const int& x, const int& y, const int& z) { assert(0); }
template<> vect<int,2>::vect (const int& x, const int& y, const int& z) { assert(0); }
template<> vect<int,4>::vect (const int& x, const int& y, const int& z) { assert(0); }

template<> vect<int,0>::vect (const int& x, const int& y, const int& z, const int& t) { assert(0); }
template<> vect<int,1>::vect (const int& x, const int& y, const int& z, const int& t) { assert(0); }
template<> vect<int,2>::vect (const int& x, const int& y, const int& z, const int& t) { assert(0); }
template<> vect<int,3>::vect (const int& x, const int& y, const int& z, const int& t) { assert(0); }



// Note: We need all dimensions all the time.
template class vect<int,0>;
template class vect<int,1>;
template class vect<int,2>;
template class vect<int,3>;
template class vect<int,4>;

template void vect<CCTK_REAL,dim>::input (istream& is);
template void vect<vect<bool,2>,dim>::input (istream& is);
template void vect<vect<bool,dim>,2>::input (istream& is);
template void vect<vect<int,dim>,2>::input (istream& is);

template void vect<bool,2>::output (ostream& os) const;
template void vect<bool,dim>::output (ostream& os) const;
template void vect<CCTK_REAL,2>::output (ostream& os) const;
template void vect<CCTK_REAL,dim>::output (ostream& os) const;
template void vect<vect<bool,2>,dim>::output (ostream& os) const;
template void vect<vect<int,2>,dim>::output (ostream& os) const;
template void vect<vect<bool,dim>,2>::output (ostream& os) const;
template void vect<vect<int,dim>,2>::output (ostream& os) const;
