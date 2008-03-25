
#include <cassert>
#include <iostream>
#include <limits>
#include <typeinfo>

#include "cctk.h"

#include "defs.hh"
#include "vect.hh"

#include "bbox.hh"

using namespace std;



// Consistency checks
template<class T, int D>
void bbox<T,D>::assert_bbox_limits () const
{
  assert (all(_stride>T(0)));
  assert (all((_upper-_lower)%_stride == T(0)));
  if (numeric_limits<T>::is_integer) {
    // prevent accidental wrap-around
    if (any (_lower >= numeric_limits<T>::max() / 2) or
        any (_lower <= numeric_limits<T>::min() / 2) or
        any (_upper >= numeric_limits<T>::max() / 2) or
        any (_upper <= numeric_limits<T>::min() / 2))
    {
      CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Tried to create a very large bbox of type %s -- it is likely that this would lead to an integer overflow",
                  typeid(*this).name());
    }
  }
}



// Accessors
template<class T, int D>
typename bbox<T,D>::size_type bbox<T,D>::size () const {
  if (empty()) return 0;
  const vect<T,D> sh(shape()/stride());
#ifdef NDEBUG
  return prod(vect<size_type,D>(sh));
#else
  size_type sz = 1, max = numeric_limits<size_type>::max();
  for (int d=0; d<D; ++d) {
    if (sh[d] > max) {
      CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "size of bbox of type %s is too large -- integer overflow",
                  typeid(*this).name());
    }
    sz *= sh[d];
    max /= sh[d];
  }
  return sz;
#endif
}

// Queries

// Containment
template<class T, int D>
bool bbox<T,D>::contains (const vect<T,D>& x) const {
  return all(x>=lower() and x<=upper());
}

template<class T, int D>
bool bbox<T,D>::is_contained_in (const bbox& b) const {
  if (empty()) return true;
  // no alignment check
  return all(lower()>=b.lower() and upper()<=b.upper());
}

// Intersection
template<class T, int D>
bool bbox<T,D>::intersects (const bbox& b) const {
  if (empty()) return false;
  if (b.empty()) return false;
  // no alignment check
  return all(upper()>=b.lower() and lower()<=b.upper());
}

// Alignment check
template<class T, int D>
bool bbox<T,D>::is_aligned_with (const bbox& b) const {
  return all(stride()==b.stride() and (lower()-b.lower()) % stride() == T(0));
}

// Operators
template<class T, int D>
bool bbox<T,D>::operator== (const bbox& b) const {
  if (empty() and b.empty()) return true;
  assert (all(stride()==b.stride()));
  return all(lower()==b.lower() and upper()==b.upper());
}

template<class T, int D>
bool bbox<T,D>::operator!= (const bbox& b) const {
  return not (*this == b);
}

#if 0
// Introduce an ordering on bboxes
template<class T, int D>
bool bbox<T,D>::operator< (const bbox& b) const {
  // An arbitraty order: empty boxes come first, then sorted by lower
  // bound, then by upper bound, then by coarseness
  if (b.empty()) return false;
  if (empty()) return true;
  for (int d=D-1; d>=0; --d) {
    if (lower()[d] < b.lower()[d]) return true;
    if (lower()[d] > b.lower()[d]) return false;
  }
  for (int d=D-1; d>=0; --d) {
    if (upper()[d] < b.upper()[d]) return true;
    if (upper()[d] > b.upper()[d]) return false;
  }
  for (int d=D-1; d>=0; --d) {
    if (stride()[d] > b.stride()[d]) return true;
    if (stride()[d] < b.stride()[d]) return false;
  }
  return false;
}
#endif

template<class T, int D>
bool bbox<T,D>::operator<= (const bbox& b) const {
  return is_contained_in (b);
}

template<class T, int D>
bool bbox<T,D>::operator>= (const bbox& b) const {
  return b <= *this;
}

template<class T, int D>
bool bbox<T,D>::operator< (const bbox& b) const {
  return *this <= b and *this != b;
}

template<class T, int D>
bool bbox<T,D>::operator> (const bbox& b) const {
  return b < *this;
}

// Expand the bbox a little by multiples of the stride
template<class T, int D>
bbox<T,D> bbox<T,D>::expand (const vect<T,D>& lo, const vect<T,D>& hi) const {
  // Allow expansion only into directions where the extent is not negative
  assert (all(lower()<=upper() or (lo==T(0) and hi==T(0))));
  const vect<T,D> str = stride();
  const vect<T,D> lb = lower() - lo * str;
  const vect<T,D> ub = upper() + hi * str;
  return bbox(lb,ub,str);
}

// Find the smallest b-compatible box around *this
template<class T, int D>
bbox<T,D> bbox<T,D>::expanded_for (const bbox& b) const {
  if (empty()) return bbox(b.lower(), b.lower()-b.stride(), b.stride());
  const vect<T,D> str = b.stride();
  const vect<T,D> loff = ((lower() - b.lower()) % str + str) % str;
  const vect<T,D> uoff = ((upper() - b.lower()) % str + str) % str;
  const vect<T,D> lo = lower() - loff; // go outwards
  const vect<T,D> up = upper() + (str - uoff) % str;
  return bbox(lo,up,str);
}

// Find the largest b-compatible box inside *this
template<class T, int D>
bbox<T,D> bbox<T,D>::contracted_for (const bbox& b) const {
  if (empty()) return bbox(b.lower(), b.lower()-b.stride(), b.stride());
  const vect<T,D> str = b.stride();
  const vect<T,D> loff = ((lower() - b.lower()) % str + str) % str;
  const vect<T,D> uoff = ((upper() - b.lower()) % str + str) % str;
  const vect<T,D> lo = lower() + (str - loff) % str; // go inwards
  const vect<T,D> up = upper() - uoff;
  return bbox(lo,up,str);
}

// Smallest bbox containing both boxes
template<class T, int D>
bbox<T,D> bbox<T,D>::expanded_containing (const bbox& b) const {
  if (empty()) return b;
  if (b.empty()) return *this;
  assert (is_aligned_with(b));
  const vect<T,D> lo = min(lower(), b.lower());
  const vect<T,D> up = max(upper(), b.upper());
  const vect<T,D> str = min(stride(), b.stride());
  return bbox(lo,up,str);
}

// Iterators
template<class T, int D>
bbox<T,D>::iterator::iterator (const bbox& box_, const vect<T,D>& pos_)
  : box(box_), pos(pos_) {
  if (box.empty()) pos=box.upper();
}

template<class T, int D>
bool bbox<T,D>::iterator::operator!= (const iterator& i) const {
  return any(pos!=i.pos);
}

template<class T, int D>
typename bbox<T,D>::iterator& bbox<T,D>::iterator::operator++ () {
  for (int d=0; d<D; ++d) {
    pos[d]+=box.stride()[d];
    if (pos[d]<=box.upper()[d]) break;
    pos[d]=box.lower()[d];
  }
  return *this;
}

template<class T, int D>
typename bbox<T,D>::iterator bbox<T,D>::begin () const {
  return iterator(*this, lower());
}

template<class T, int D>
typename bbox<T,D>::iterator bbox<T,D>::end () const {
  return iterator(*this, lower());
}



// Input
template<class T,int D>
void bbox<T,D>::input (istream& is) {
  try {
    skipws (is);
    consume (is, '(');
    is >> _lower;
    skipws (is);
    consume (is, ':');
    is >> _upper;
    skipws (is);
    consume (is, ':');
    is >> _stride;
    skipws (is);
    if (is.peek() == '/') {
      consume (is, '/');
      vect<T,D> lower_dummy;
      is >> lower_dummy;
      skipws (is);
      consume (is, ':');
      vect<T,D> upper_dummy;
      is >> upper_dummy;
      skipws (is);
      consume (is, '/');
      vect<T,D> shape_dummy;
      is >> shape_dummy;
      skipws (is);
      consume (is, '/');
      size_type size_dummy;
      is >> size_dummy;
      assert (is.good());
      skipws (is);
    }
    consume (is, ')');
  } catch (input_error &err) {
    cout << "Input error while reading a bbox" << endl;
    throw err;
  }
  if (any(_stride<=T(0))) {
    cout << "While reading the bbox " << *this << ":" << endl
         << "   The stride is not positive." << endl;
    throw input_error();
  }
  if (any((_upper-_lower)%_stride != T(0))) {
    cout << "While reading the bbox " << *this << ":" << endl
         << "   The stride does not evenly divide the extent." << endl;
    throw input_error();
  }
  assert (all(_stride>T(0)));
  assert (all((_upper-_lower)%_stride == T(0)));
}



// Output
template<class T,int D>
void bbox<T,D>::output (ostream& os) const {
  os << "(" << lower() << ":" << upper() << ":" << stride()
     << "/" << lower() / stride() << ":" << upper() / stride()
     << "/" << shape() / stride()
     << "/" << size() << ")";
}



// Note: We need all dimensions all the time.
template class bbox<int,0>;
template class bbox<int,1>;
template class bbox<int,2>;
template class bbox<int,3>;
template class bbox<CCTK_REAL,3>;
