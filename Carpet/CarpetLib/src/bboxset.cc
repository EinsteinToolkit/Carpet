// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/bboxset.cc,v 1.12 2003/02/25 22:57:00 schnetter Exp $

#include <assert.h>

#include <iostream>
#include <set>

#include "defs.hh"

#include "bboxset.hh"

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

template<class T, int D>
bboxset<T,D> bboxset<T,D>::plus (const bbox<T,D>& b1, const bbox<T,D>& b2) {
  return bboxset(b1) + b2;
}

template<class T, int D>
bboxset<T,D> bboxset<T,D>::plus (const bbox<T,D>& b, const bboxset<T,D>& s) {
  return s + b;
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
template<class T, int D>
bboxset<T,D> bboxset<T,D>::minus (const bbox<T,D>& b1, const bbox<T,D>& b2) {
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

template<class T, int D>
bboxset<T,D> bboxset<T,D>::minus (const bbox<T,D>& b, const bboxset<T,D>& s) {
  bboxset<T,D> r = bboxset<T,D>(b) - s;
  assert (r.invariant());
  return r;
}



// Equality
template<class T, int D>
bool bboxset<T,D>::operator<= (const bboxset<T,D>& s) const {
  return (*this - s).empty();
}

template<class T, int D>
bool bboxset<T,D>::operator< (const bboxset<T,D>& s) const {
  return (*this - s).empty() && ! (s - *this).empty();
}

template<class T, int D>
bool bboxset<T,D>::operator>= (const bboxset<T,D>& s) const {
  return s <= *this;
}

template<class T, int D>
bool bboxset<T,D>::operator> (const bboxset<T,D>& s) const {
  return s < *this;
}

template<class T, int D>
bool bboxset<T,D>::operator== (const bboxset<T,D>& s) const {
  return (*this <= s) && (*this >= s);
}

template<class T, int D>
bool bboxset<T,D>::operator!= (const bboxset<T,D>& s) const {
  return ! (*this == s);
}



// Output
template<class T,int D>
void bboxset<T,D>::output (ostream& os) const {
  T Tdummy;
  os << "bboxset<" << typestring(Tdummy) << "," << D << ">:"
     << "size=" << size() << "," << "set=" << bs;
}



template class bboxset<int,3>;
