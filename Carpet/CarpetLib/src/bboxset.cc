/***************************************************************************
                          bboxset.cc  -  Sets of bounding boxes
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/bboxset.cc,v 1.5 2001/03/22 18:42:05 eschnett Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <assert.h>

#include <iostream>
#include <set>

#include "defs.hh"

#if !defined(TMPL_IMPLICIT) && !defined(BBOXSET_HH)
#  include "bboxset.hh"
#endif

using namespace std;



// Constructors
template<class T, int D>
bboxset<T,D>::bboxset () {
  assert (invariant());
}

template<class T, int D>
bboxset<T,D>::bboxset (const box& b) {
  if (!b.empty()) bs.insert(b);
  assert (invariant());
}

template<class T, int D>
bboxset<T,D>::bboxset (const bboxset& s): bs(s.bs) {
  assert (invariant());
}

template<class T, int D>
bboxset<T,D>::bboxset (const bset& bs): bs(bs) {
  assert (invariant());
}



// Invariant
template<class T, int D>
bool bboxset<T,D>::invariant () const {
  for (const_iterator bi=begin(); bi!=end(); ++bi) {
    if ((*bi).empty()) return false;
    if (! (*bi).is_aligned_with(*bs.begin())) return false;
    // check for overlap (quadratic -- expensive)
    int cnt=0;
    for (const_iterator bi2=bi; bi2!=end(); ++bi2) {
      if (!cnt++) continue;
      if (! ((*bi2) & (*bi)).empty()) return false;
    }
  }
  return true;
}



// Normalisation
template<class T, int D>
void bboxset<T,D>::normalize () {
  // TODO
  assert (invariant());
}



// Accessors
template<class T, int D>
T bboxset<T,D>::size () const {
  T s=0;
  for (const_iterator bi=begin(); bi!=end(); ++bi) {
    s += (*bi).size();
  }
  return s;
}



// Add (bboxes that don't overlap)
template<class T, int D>
bboxset<T,D>& bboxset<T,D>::operator+= (const box& b) {
  if (b.empty()) return *this;
  // check for overlap
  for (const_iterator bi=begin(); bi!=end(); ++bi) {
    assert ((*bi & b).empty());
  }
  bs.insert(b);
  assert (invariant());
  return *this;
}

template<class T, int D>
bboxset<T,D>& bboxset<T,D>::operator+= (const bboxset& s) {
  for (const_iterator bi=s.begin(); bi!=s.end(); ++bi) {
    *this += *bi;
  }
  assert (invariant());
  return *this;
}

template<class T, int D>
bboxset<T,D> bboxset<T,D>::operator+ (const box& b) const {
  bboxset r(*this);
  r += b;
  assert (r.invariant());
  return r;
}

template<class T, int D>
bboxset<T,D> bboxset<T,D>::operator+ (const bboxset& s) const {
  bboxset r(*this);
  r += s;
  assert (r.invariant());
  return r;
}



// Union
template<class T, int D>
bboxset<T,D>& bboxset<T,D>::operator|= (const box& b) {
  *this += b - *this;
  assert (invariant());
  return *this;
}

template<class T, int D>
bboxset<T,D>& bboxset<T,D>::operator|= (const bboxset& s) {
  *this += s - *this;
  assert (invariant());
  return *this;
}

template<class T, int D>
bboxset<T,D> bboxset<T,D>::operator| (const box& b) const {
  bboxset r(*this);
  r |= b;
  assert (r.invariant());
  return r;
}

template<class T, int D>
bboxset<T,D> bboxset<T,D>::operator| (const bboxset& s) const {
  bboxset r(*this);
  r |= s;
  assert (r.invariant());
  return r;
}



// Intersection
template<class T, int D>
bboxset<T,D> bboxset<T,D>::operator& (const box& b) const {
  // start with an empty set
  bboxset r;
  // walk all my elements
  for (const_iterator bi=begin(); bi!=end(); ++bi) {
    // insert the intersection with the bbox
    r += *bi & b;
  }
  assert (r.invariant());
  return r;
}

template<class T, int D>
bboxset<T,D> bboxset<T,D>::operator& (const bboxset& s) const {
  // start with an empty set
  bboxset r;
  // walk all the bboxes
  for (const_iterator bi=s.begin(); bi!=s.end(); ++bi) {
    // insert the intersection with this bbox
    r += *this & *bi;
  }
  assert (r.invariant());
  return r;
}

template<class T, int D>
bboxset<T,D>& bboxset<T,D>::operator&= (const box& b) {
  *this = *this & b;
  assert (invariant());
  return *this;
}

template<class T, int D>
bboxset<T,D>& bboxset<T,D>::operator&= (const bboxset& s) {
  *this = *this & s;
  assert (invariant());
  return *this;
}



// Difference
#ifndef SGI
// This doesn't work on SGIs.  Is this legal C++?
template<class T, int D>
bboxset<T,D> operator- (const bbox<T,D>& b1, const bbox<T,D>& b2) {
  assert (b1.is_aligned_with(b2));
  if (b1.empty()) return bboxset<T,D>();
  if (b2.empty()) return bboxset<T,D>(b1);
  const vect<T,D> str = b1.stride();
  bboxset<T,D> r;
  for (int d=0; d<D; ++d) {
    // make resulting bboxes as large as possible in x-direction (for
    // better consumption by Fortranly ordered arrays)
    vect<T,D> lb, ub;
    bbox<T,D> b;
    for (int dd=0; dd<D; ++dd) {
      if (dd<d) {
	lb[dd] = b2.lower()[dd];
	ub[dd] = b2.upper()[dd];
      } else if (dd>d) {
	lb[dd] = b1.lower()[dd];
	ub[dd] = b1.upper()[dd];
      }
    }
    lb[d] = b1.lower()[d];
    ub[d] = b2.lower()[d] - str[d];
    b = bbox<T,D>(lb,ub,str) & b1;
    r += b;
    lb[d] = b2.upper()[d] + str[d];
    ub[d] = b1.upper()[d];
    b = bbox<T,D>(lb,ub,str) & b1;
    r += b;
  }
  assert (r.invariant());
  return r;
}
#else
bboxset<int,1> operator- (const bbox<int,1>& b1, const bbox<int,1>& b2) {
  assert (b1.is_aligned_with(b2));
  if (b1.empty()) return bboxset<int,1>();
  if (b2.empty()) return bboxset<int,1>(b1);
  const vect<int,1> str = b1.stride();
  bboxset<int,1> r;
  for (int d=0; d<1; ++d) {
    // make resulting bboxes as large as possible in x-direction (for
    // better consumption by Fortranly ordered arrays)
    vect<int,1> lb, ub;
    bbox<int,1> b;
    for (int dd=0; dd<1; ++dd) {
      if (dd<d) {
	lb[dd] = b2.lower()[dd];
	ub[dd] = b2.upper()[dd];
      } else if (dd>d) {
	lb[dd] = b1.lower()[dd];
	ub[dd] = b1.upper()[dd];
      }
    }
    lb[d] = b1.lower()[d];
    ub[d] = b2.lower()[d] - str[d];
    b = bbox<int,1>(lb,ub,str) & b1;
    r += b;
    lb[d] = b2.upper()[d] + str[d];
    ub[d] = b1.upper()[d];
    b = bbox<int,1>(lb,ub,str) & b1;
    r += b;
  }
  assert (r.invariant());
  return r;
}
bboxset<int,2> operator- (const bbox<int,2>& b1, const bbox<int,2>& b2) {
  assert (b1.is_aligned_with(b2));
  if (b1.empty()) return bboxset<int,2>();
  if (b2.empty()) return bboxset<int,2>(b1);
  const vect<int,2> str = b1.stride();
  bboxset<int,2> r;
  for (int d=0; d<2; ++d) {
    // make resulting bboxes as large as possible in x-direction (for
    // better consumption by Fortranly ordered arrays)
    vect<int,2> lb, ub;
    bbox<int,2> b;
    for (int dd=0; dd<2; ++dd) {
      if (dd<d) {
	lb[dd] = b2.lower()[dd];
	ub[dd] = b2.upper()[dd];
      } else if (dd>d) {
	lb[dd] = b1.lower()[dd];
	ub[dd] = b1.upper()[dd];
      }
    }
    lb[d] = b1.lower()[d];
    ub[d] = b2.lower()[d] - str[d];
    b = bbox<int,2>(lb,ub,str) & b1;
    r += b;
    lb[d] = b2.upper()[d] + str[d];
    ub[d] = b1.upper()[d];
    b = bbox<int,2>(lb,ub,str) & b1;
    r += b;
  }
  assert (r.invariant());
  return r;
}
bboxset<int,3> operator- (const bbox<int,3>& b1, const bbox<int,3>& b2) {
  assert (b1.is_aligned_with(b2));
  if (b1.empty()) return bboxset<int,3>();
  if (b2.empty()) return bboxset<int,3>(b1);
  const vect<int,3> str = b1.stride();
  bboxset<int,3> r;
  for (int d=0; d<3; ++d) {
    // make resulting bboxes as large as possible in x-direction (for
    // better consumption by Fortranly ordered arrays)
    vect<int,3> lb, ub;
    bbox<int,3> b;
    for (int dd=0; dd<3; ++dd) {
      if (dd<d) {
	lb[dd] = b2.lower()[dd];
	ub[dd] = b2.upper()[dd];
      } else if (dd>d) {
	lb[dd] = b1.lower()[dd];
	ub[dd] = b1.upper()[dd];
      }
    }
    lb[d] = b1.lower()[d];
    ub[d] = b2.lower()[d] - str[d];
    b = bbox<int,3>(lb,ub,str) & b1;
    r += b;
    lb[d] = b2.upper()[d] + str[d];
    ub[d] = b1.upper()[d];
    b = bbox<int,3>(lb,ub,str) & b1;
    r += b;
  }
  assert (r.invariant());
  return r;
}
#endif

template<class T, int D>
bboxset<T,D> bboxset<T,D>::operator- (const box& b) const {
  // start with an empty set
  bboxset r;
  // walk all my elements
  for (const_iterator bi=begin(); bi!=end(); ++bi) {
    // insert the difference with the bbox
    r += *bi - b;
  }
  assert (r.invariant());
  return r;
}

template<class T, int D>
bboxset<T,D>& bboxset<T,D>::operator-= (const box& b) {
  *this = *this - b;
  assert (invariant());
  return *this;
}
  
template<class T, int D>
bboxset<T,D>& bboxset<T,D>::operator-= (const bboxset& s) {
  for (const_iterator bi=s.begin(); bi!=s.end(); ++bi) {
    *this -= *bi;
  }
  assert (invariant());
  return *this;
}

template<class T, int D>
bboxset<T,D> bboxset<T,D>::operator- (const bboxset& s) const {
  bboxset r(*this);
  r -= s;
  assert (r.invariant());
  return r;
}

#ifndef SGI
// This doesn't work on SGIs.  Is this legal C++?
template<class T, int D>
bboxset<T,D> operator- (const bbox<T,D>& b, const bboxset<T,D>& s) {
  bboxset<T,D> r = bboxset<T,D>(b) - s;
  assert (r.invariant());
  return r;
}
#else
bboxset<int,1> operator- (const bbox<int,1>& b, const bboxset<int,1>& s) {
  bboxset<int,1> r = bboxset<int,1>(b) - s;
  assert (r.invariant());
  return r;
}
bboxset<int,2> operator- (const bbox<int,2>& b, const bboxset<int,2>& s) {
  bboxset<int,2> r = bboxset<int,2>(b) - s;
  assert (r.invariant());
  return r;
}
bboxset<int,3> operator- (const bbox<int,3>& b, const bboxset<int,3>& s) {
  bboxset<int,3> r = bboxset<int,3>(b) - s;
  assert (r.invariant());
  return r;
}
#endif



// Output
template<class T,int D>
ostream& operator<< (ostream& os, const bboxset<T,D>& s) {
  os << "bboxset<T," << D << ">:size=" << s.size() << "," << "set=" << s.bs;
  return os;
}



#if defined(TMPL_EXPLICIT)
template class bboxset<int,3>;

template
bboxset<int,3> operator- (const bbox<int,3>& b1, const bbox<int,3>& b3);
template
bboxset<int,3> operator- (const bbox<int,3>& b, const bboxset<int,3>& s);
template
ostream& operator<< (ostream& os, const bboxset<int,3>& b);
#endif
