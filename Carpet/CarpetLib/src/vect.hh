/***************************************************************************
                          vect.hh  -  Small inline vectors
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/vect.hh,v 1.2 2001/03/17 22:37:56 eschnett Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef VECT_HH
#define VECT_HH

#include <cassert>
#include <cmath>
#include <iostream>



// Forward definition
template<class T, int D> class vect;

// Output
template<class T,int D>
ostream& operator<< (ostream& os, const vect<T,D>& a);



// Vect class
template<class T, int D>
class vect {
  
  // Fields
  T elt[D];
  
public:
  
  // Constructors
  explicit vect () { }
  
  vect (const vect& a) {
    for (int d=0; d<D; ++d) elt[d]=a.elt[d];
  }
  
  // this constructor might be confusing, but is so convenient
  vect (const T x) {
    for (int d=0; d<D; ++d) elt[d]=x;
  }
  
  vect (const T x, const T y) {
    assert (D==2);
    elt[0]=x; elt[1]=y;
  }
  
  vect (const T x, const T y, const T z) {
    assert (D==3);
    elt[0]=x; elt[1]=y; elt[2]=z;
  }
  
  vect (const T* const x) {
    for (int d=0; d<D; ++d) elt[d]=x[d];
  }
  
  template<class S>
  explicit vect (const vect<S,D>& a) {
    for (int d=0; d<D; ++d) elt[d]=(T)a[d];
  }
  
  static vect dir (const int d) {
    vect r=0;
    r[d]=1;
    return r;
  }
  
  static vect seq () {
    vect r;
    for (int d=0; d<D; ++d) r[d]=d;
    return r;
  }
  
  static vect seq (const int n) {
    vect r;
    for (int d=0; d<D; ++d) r[d]=n+d;
    return r;
  }
  
  static vect seq (const int n, const int s) {
    vect r;
    for (int d=0; d<D; ++d) r[d]=n+s*d;
    return r;
  }
  
  // Accessors
  const T& operator[] (const int d) const {
    assert(d>=0 && d<D);
    return elt[d];
  }
  
  T& operator[] (const int d) {
    assert(d>=0 && d<D);
    return elt[d];
  }
  
  template<class TT, int DD>
  vect<T,DD> operator[] (const vect<TT,DD>& a) const {
    vect<T,DD> r;
    for (int d=0; d<DD; ++d) r[d] = (*this)[a[d]];
    return r;
  }
  
  // Modifying operators
  vect& operator+=(const T x) {
    for (int d=0; d<D; ++d) (*this)[d]+=x;
    return *this;
  }
  
  vect& operator-=(const T x) {
    for (int d=0; d<D; ++d) (*this)[d]-=x;
    return *this;
  }
  
  vect& operator*=(const T x) {
    for (int d=0; d<D; ++d) (*this)[d]*=x;
    return *this;
  }
  
  vect& operator/=(const T x) {
    for (int d=0; d<D; ++d) (*this)[d]/=x;
    return *this;
  }
  
  vect& operator%=(const T x) {
    for (int d=0; d<D; ++d) (*this)[d]%=x;
    return *this;
  }
  
  vect& operator&=(const T x) {
    for (int d=0; d<D; ++d) (*this)[d]&=x;
    return *this;
  }
  
  vect& operator|=(const T x) {
    for (int d=0; d<D; ++d) (*this)[d]|=x;
    return *this;
  }
  
  vect& operator^=(const T x) {
    for (int d=0; d<D; ++d) (*this)[d]^=x;
    return *this;
  }
  
  vect& operator+=(const vect& a) {
    for (int d=0; d<D; ++d) (*this)[d]+=a[d];
    return *this;
  }
  
  vect& operator-=(const vect& a) {
    for (int d=0; d<D; ++d) (*this)[d]-=a[d];
    return *this;
  }
  
  vect& operator*=(const vect& a) {
    for (int d=0; d<D; ++d) (*this)[d]*=a[d];
    return *this;
  }
  
  vect& operator/=(const vect& a) {
    for (int d=0; d<D; ++d) (*this)[d]/=a[d];
    return *this;
  }
  
  vect& operator%=(const vect& a) {
    for (int d=0; d<D; ++d) (*this)[d]%=a[d];
    return *this;
  }
  
  vect& operator&=(const vect& a) {
    for (int d=0; d<D; ++d) (*this)[d]&=a[d];
    return *this;
  }
  
  vect& operator|=(const vect& a) {
    for (int d=0; d<D; ++d) (*this)[d]|=a[d];
    return *this;
  }
  
  vect& operator^=(const vect& a) {
    for (int d=0; d<D; ++d) (*this)[d]^=a[d];
    return *this;
  }
  
  // Non-modifying operators
  vect operator+ () const {
    vect r(*this);
    for (int d=0; d<D; ++d) r[d]=+r[d];
    return r;
  }
  
  vect operator- () const {
    vect r(*this);
    for (int d=0; d<D; ++d) r[d]=-r[d];
    return r;
  }
  
  vect<bool,D> operator! () const {
    vect r(*this);
    for (int d=0; d<D; ++d) r[d]=!r[d];
    return r;
  }
  
  vect operator~ () const {
    vect r(*this);
    for (int d=0; d<D; ++d) r[d]=~r[d];
    return r;
  }
  
  vect operator+ (const T x) const {
    vect r(*this);
    r+=x;
    return r;
  }
  
  vect operator- (const T x) const {
    vect r(*this);
    r-=x;
    return r;
  }
  
  vect operator* (const T x) const {
    vect r(*this);
    r*=x;
    return r;
  }
  
  vect operator/ (const T x) const {
    vect r(*this);
    r/=x;
    return r;
  }
  
  vect operator% (const T x) const {
    vect r(*this);
    r%=x;
    return r;
  }
  
  vect operator& (const T x) const {
    vect r(*this);
    r&=x;
    return r;
  }
  
  vect operator| (const T x) const {
    vect r(*this);
    r|=x;
    return r;
  }
  
  vect operator^ (const T x) const {
    vect r(*this);
    r^=x;
    return r;
  }
  
  vect<bool,D> operator&& (const T x) const {
    vect r(*this);
    for (int d=0; d<D; ++d) r[d]=r[d]&&x;
    return r;
  }
  
  vect<bool,D> operator|| (const T x) const {
    vect r(*this);
    for (int d=0; d<D; ++d) r[d]=r[d]||x;
    return r;
  }
  
  vect operator+ (const vect& a) const {
    vect r(*this);
    r+=a;
    return r;
  }
  
  vect operator- (const vect& a) const {
    vect r(*this);
    r-=a;
    return r;
  }
  
  vect operator* (const vect& a) const {
    vect r(*this);
    r*=a;
    return r;
  }
  
  vect operator/ (const vect& a) const {
    vect r(*this);
    r/=a;
    return r;
  }
  
  vect operator% (const vect& a) const {
    vect r(*this);
    r%=a;
    return r;
  }
  
  vect operator& (const vect& a) const {
    vect r(*this);
    r&=a;
    return r;
  }
  
  vect operator| (const vect& a) const {
    vect r(*this);
    r|=a;
    return r;
  }
  
  vect operator^ (const vect& a) const {
    vect r(*this);
    r^=a;
    return r;
  }
  
  vect<bool,D> operator&& (const vect& a) const {
    vect<bool,D> r;
    for (int d=0; d<D; ++d) r[d]=(*this)[d]&&a[d];
    return r;
  }
  
  vect<bool,D> operator|| (const vect& a) const {
    vect<bool,D> r;
    for (int d=0; d<D; ++d) r[d]=(*this)[d]||a[d];
    return r;
  }
  
  vect<bool,D> operator== (const vect& a) const {
    vect<bool,D> r;
    for (int d=0; d<D; ++d) r[d]=(*this)[d]==a[d];
    return r;
  }
  
  vect<bool,D> operator!= (const vect& a) const {
    vect<bool,D> r;
    for (int d=0; d<D; ++d) r[d]=(*this)[d]!=a[d];
    return r;
  }
  
  vect<bool,D> operator< (const vect& a) const {
    vect<bool,D> r;
    for (int d=0; d<D; ++d) r[d]=(*this)[d]<a[d];
    return r;
  }
  
  vect<bool,D> operator<= (const vect& a) const {
    vect<bool,D> r;
    for (int d=0; d<D; ++d) r[d]=(*this)[d]<=a[d];
    return r;
  }
  
  vect<bool,D> operator> (const vect& a) const {
    vect<bool,D> r;
    for (int d=0; d<D; ++d) r[d]=(*this)[d]>a[d];
    return r;
  }
  
  vect<bool,D> operator>= (const vect& a) const {
    vect<bool,D> r;
    for (int d=0; d<D; ++d) r[d]=(*this)[d]>=a[d];
    return r;
  }
  
#if 0
  vect operator?: (const vect& a, const vect& b) const {
    vect r(*this);
    for (int d=0; d<D; ++d) r[d]=r[d]?a[d]:b[d];
    return r;
  }
#endif
  
  // Iterators
#if 0
  // This is non-standard
  class iter {
    vect &vec;
    int d;
  public:
    iter (vect &a): vec(a), d(0) { }
    iter& operator++ () { assert(d<D); ++d; return *this; }
    bool operator bool () { return d==D; }
    T& operator* { return vec[d]; }
  };
#endif
  
  // Output
  friend ostream& operator<< <>(ostream& os, const vect& a);
};



// Operators
template<class T,int D>
inline vect<T,D> abs (const vect<T,D>& a) {
  vect<T,D> r;
  for (int d=0; d<D; ++d) r[d]=abs(a[d]);
  return r;
}

template<class T,int D>
inline vect<T,D> max (const vect<T,D>& a, const vect<T,D>& b) {
  vect<T,D> r;
  for (int d=0; d<D; ++d) r[d]=max(a[d],b[d]);
  return r;
}

template<class T,int D>
inline vect<T,D> min (const vect<T,D>& a, const vect<T,D>& b) {
  vect<T,D> r;
  for (int d=0; d<D; ++d) r[d]=min(a[d],b[d]);
  return r;
}



// Reduction operators
template<int D>
inline bool any (const vect<bool,D>& a) {
  bool r(false);
  for (int d=0; d<D; ++d) r|=a[d];
  return r;
}

template<int D>
inline bool all (const vect<bool,D>& a) {
  bool r(true);
  for (int d=0; d<D; ++d) r&=a[d];
  return r;
}

template<class T,int D>
inline int count (const vect<T,D>& a) {
  return D;
}

template<class T,int D>
inline T dot (const vect<T,D>& a, const vect<T,D>& b) {
  T r(0);
  for (int d=0; d<D; ++d) r+=a[d]*b[d];
  return r;
}

template<class T,int D>
inline T hypot (const vect<T,D>& a) {
  return sqrt(dot(a,a));
}

template<class T,int D>
inline T maxval (const vect<T,D>& a) {
  assert (D>0);
  T r(a[0]);
  for (int d=0; d<D; ++d) r=max(r,a[d]);
  return r;
}

template<class T,int D>
inline T minval (const vect<T,D>& a) {
  assert (D>0);
  T r(a[0]);
  for (int d=0; d<D; ++d) r=min(r,a[d]);
  return r;
}

template<class T,int D>
inline T prod (const vect<T,D>& a) {
  T r(1);
  for (int d=0; d<D; ++d) r*=a[d];
  return r;
}

template<class T,int D>
inline int size (const vect<T,D>& a) {
  return D;
}

template<class T,int D>
inline T sum (const vect<T,D>& a) {
  T r(0);
  for (int d=0; d<D; ++d) r+=a[d];
  return r;
}



#if defined(TMPL_IMPLICIT)
#  include "vect.cc"
#endif

#endif // VECT_HH
