/***************************************************************************
                          bbox.cc  -  Bounding boxes
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/bbox.cc,v 1.3 2001/03/07 13:00:57 eschnett Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <cassert>
#include <iostream>

#include "defs.hh"
#include "vect.hh"

#if !defined(TMPL_IMPLICIT) || !defined(BBOX_HH)
#  include "bbox.hh"
#endif



// Constructors
template<class T, int D>
bbox<T,D>::bbox (): _lower(1), _upper(0), _stride(1) { }

template<class T, int D>
bbox<T,D>::bbox (const bbox& b)
  : _lower(b._lower), _upper(b._upper), _stride(b._stride)
{ }

template<class T, int D>
bbox<T,D>& bbox<T,D>::operator= (const bbox& b) {
  _lower=b._lower; _upper=b._upper; _stride=b._stride;
  return *this;
}

template<class T, int D>
bbox<T,D>::bbox (const vect<T,D>& lower, const vect<T,D>& upper,
		 const vect<T,D>& stride)
  : _lower(lower), _upper(upper), _stride(stride)
{
  assert (all(stride>=1));
  assert (all((upper-lower)%stride==0));
}

// Accessors
template<class T, int D>
T bbox<T,D>::size () const {
  if (empty()) return 0;
  return prod(shape());
}

template<class T, int D>
T bbox<T,D>::num_points () const {
  if (empty()) return 0;
  return prod((shape()+stride()-1)/stride());
}

// Queries
template<class T, int D>
bool bbox<T,D>::contains (const vect<T,D>& x) const {
  return all(x>=lower() && x<=upper());
}

// Operators
template<class T, int D>
bool bbox<T,D>::operator== (const bbox& b) const {
  if (empty() && b.empty()) return true;
  assert (all(stride()==b.stride()));
  return all(lower()==b.lower() && upper()==b.upper());
}

template<class T, int D>
bool bbox<T,D>::operator!= (const bbox& b) const {
  return ! (*this == b);
}

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

template<class T, int D>
bool bbox<T,D>::operator> (const bbox& b) const {
  return b < *this;
}

template<class T, int D>
bool bbox<T,D>::operator<= (const bbox& b) const {
  return ! (b > *this);
}

template<class T, int D>
bool bbox<T,D>::operator>= (const bbox& b) const {
  return b <= *this;
}

// Intersection
template<class T, int D>
bbox<T,D> bbox<T,D>::operator& (const bbox& b) const {
  assert (all(stride()==b.stride()));
  vect<T,D> lo = max(lower(),b.lower());
  vect<T,D> up = min(upper(),b.upper());
  return bbox(lo,up,stride());
}

// Containment
template<class T, int D>
bool bbox<T,D>::contained_in (const bbox& b) const {
  // no alignment check
  return all(lower()>=b.lower() && upper()<=b.upper());
}

// Alignment check
template<class T, int D>
bool bbox<T,D>::aligned_with (const bbox& b) const {
  return all(stride()==b.stride() && (lower()-b.lower()) % stride() == 0);
}

// Expand the bbox a little by multiples of the stride
template<class T, int D>
bbox<T,D> bbox<T,D>::expand (const vect<T,D>& lo, const vect<T,D>& hi) const {
  const vect<T,D> str = stride();
  const vect<T,D> lb = lower() - lo * str;
  const vect<T,D> ub = upper() + hi * str;
  return bbox(lb,ub,str);
}

// Find the smallest b-compatible box around *this
template<class T, int D>
bbox<T,D> bbox<T,D>::expanded_for (const bbox& b) const {
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
  const vect<T,D> str = b.stride();
  const vect<T,D> loff = ((lower() - b.lower()) % str + str) % str;
  const vect<T,D> uoff = ((upper() - b.lower()) % str + str) % str;
  const vect<T,D> lo = lower() + (str - loff) % str; // go inwards
  const vect<T,D> up = upper() - uoff;
  return bbox(lo,up,str);
}

// Set operations
// Smallest bbox containing both boxes
template<class T, int D>
bbox<T,D> bbox<T,D>::operator* (const bbox& b) const {
  if (empty()) return b;
  if (b.empty()) return *this;
  assert (aligned_with(b));
  const vect<T,D> lo = min(lower(), b.lower());
  const vect<T,D> up = max(upper(), b.upper());
  const vect<T,D> str = min(stride(), b.stride());
  return bbox(lo,up,str);
}

template<class T, int D>
bbox<T,D>& bbox<T,D>::operator*= (const bbox& b) {
  *this = *this * b;
  return *this;
}

// Largest bbox inside both boxes
template<class T, int D>
bbox<T,D> bbox<T,D>::operator+ (const bbox& b) const {
  if (empty() || b.empty()) return bbox();
  assert (aligned_with(b));
  const vect<T,D> lo = max(lower(), b.lower());
  const vect<T,D> up = min(upper(), b.upper());
  const vect<T,D> str = min(stride(), b.stride());
  return bbox(lo,up,str);
}

template<class T, int D>
bbox<T,D>& bbox<T,D>::operator+= (const bbox& b) {
  *this = *this + b;
  return *this;
}

// Iterators
template<class T, int D>
bbox<T,D>::iterator::iterator (const bbox& box, const vect<T,D>& pos)
  : box(box), pos(pos) {
  if (box.empty()) this->pos=box.upper()+box.stride();
}

template<class T, int D>
bool bbox<T,D>::iterator::operator!= (const iterator& i) const {
  return any(pos!=i.pos);
}

template<class T, int D>
bbox<T,D>::iterator& bbox<T,D>::iterator::operator++ () {
  for (int d=0; d<D; ++d) {
    pos[d]+=box.stride()[d];
    if (pos[d]<=box.upper()[d]) return *this;
    pos[d]=box.lower()[d];
  }
  pos=box.end().pos;
  return *this;
}

template<class T, int D>
bbox<T,D>::iterator bbox<T,D>::begin () const {
  return iterator(*this, lower());
}

template<class T, int D>
bbox<T,D>::iterator bbox<T,D>::end () const {
  return iterator(*this, upper()+stride());
}



// Output
template<class T,int D>
ostream& operator<< (ostream& os, const bbox<T,D>& b) {
  os << "(" << b.lower() << ":" << b.upper() << ":" << b.stride() << ")";
  return os;
}



#if defined(TMPL_EXPLICIT)
template class bbox<int,1>;
template ostream& operator<< (ostream& os, const bbox<int,1>& b);
template class bbox<int,2>;
template ostream& operator<< (ostream& os, const bbox<int,2>& b);
template class bbox<int,3>;
template ostream& operator<< (ostream& os, const bbox<int,3>& b);
#endif
