#include <cassert>
#include <iostream>

#include "cctk.h"

#include "defs.hh"
#include "bboxset.hh"

#include "vect.hh"

using namespace std;



// Input
template<typename T,int D>
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
template<typename T,int D>
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
template void vect<vect<CCTK_REAL,dim>,2>::output (ostream& os) const;



// Instantiate for bboxset class

#define DEFINE_FAKE_VECT_OPERATIONS(T,D)                                \
template<> vect<T,D> vect<T,D>::dir (const int d) { assert(0); }        \
template<> vect<T,D> vect<T,D>::seq () { assert(0); }                   \
template<> vect<T,D> vect<T,D>::seq (const int n) { assert(0); }        \
template<> vect<T,D> vect<T,D>::seq (const int n, const int s) { assert(0); } \
template<> vect<T,D>& vect<T,D>::operator*= (const vect<T,D>&) { assert(0); } \
template<> vect<T,D>& vect<T,D>::operator*= (const T&) { assert(0); }   \
template<> vect<T,D>& vect<T,D>::operator/= (const vect<T,D>&) { assert(0); } \
template<> vect<T,D>& vect<T,D>::operator/= (const T&) { assert(0); }   \
template<> vect<T,D>& vect<T,D>::operator%= (const vect<T,D>&) { assert(0); } \
template<> vect<T,D>& vect<T,D>::operator%= (const T&) { assert(0); }   \
template<> vect<T,D>& vect<T,D>::operator^= (const vect<T,D>&) { assert(0); } \
template<> vect<T,D>& vect<T,D>::operator^= (const T&) { assert(0); }   \
template<> vect<T,D> vect<T,D>::operator+ () const { assert(0); }       \
template<> vect<T,D> vect<T,D>::operator- () const { assert(0); }       \
template<> vect<T,D> vect<T,D>::operator~ () const { assert(0); }       \
template class vect<T,D>;                                               \
template size_t memoryof (const vect<T,D>&);                            \
template istream& operator>> (istream& is, vect<T,D>&);                 \
template ostream& operator<< (ostream& os, const vect<T,D>&);

typedef bboxset<int,dim> T1;
typedef vect<bboxset<int,dim>,2> T2;

DEFINE_FAKE_VECT_OPERATIONS(T1,dim)
DEFINE_FAKE_VECT_OPERATIONS(T2,dim)

#undef DEFINE_FAKE_VECT_OPERATIONS
