// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/vect.hh,v 1.14 2003/03/18 17:30:26 schnetter Exp $

#ifndef VECT_HH
#define VECT_HH

#include <assert.h>
#include <math.h>

#include <iostream>

using namespace std;



// Forward definition
template<class T, int D> class vect;

// Input/Output
template<class T,int D>
istream& operator>> (istream& is, vect<T,D>& a);
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
  
  // this constructor might be confusing, but it is so convenient
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
  
  vect (const T x, const T y, const T z, const T t) {
    assert (D==4);
    elt[0]=x; elt[1]=y; elt[2]=z; elt[3]=t;
  }
  
  vect (const T* const x) {
    for (int d=0; d<D; ++d) elt[d]=x[d];
  }
  
  template<class S>
  explicit vect (const vect<S,D>& a) {
    for (int d=0; d<D; ++d) elt[d]=(T)a[d];
  }
  
  static vect make () {
    assert (D==0);
    return vect();
  }
  
  static vect make (const T x) {
    assert (D==1);
    return vect(x);
  }
  
  static vect make (const T x, const T y) {
    assert (D==2);
    return vect(x, y);
  }
  
  static vect make (const T x, const T y, const T z) {
    assert (D==3);
    return vect(x, y, z);
  }
  
  static vect make (const T x, const T y, const T z, const T t) {
    assert (D==4);
    return vect(x, y, z, t);
  }
  
  static vect& ref (T* const x) {
    return *(vect*)x;
  }
  
  static vect dir (const int d) {
    vect r=(T)0;
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
  vect replace (const int d, const T x) {
    assert (d>=0 && d<D);
    vect r;
    for (int dd=0; dd<D; ++dd) r[dd]=dd==d?x:(*this)[dd];
    return r;
  }
  
  vect operator+ () const {
    vect r;
    for (int d=0; d<D; ++d) r[d]=+r[d];
    return r;
  }
  
  vect operator- () const {
    vect r;
    for (int d=0; d<D; ++d) r[d]=-r[d];
    return r;
  }
  
  vect<bool,D> operator! () const {
    vect<bool,D> r;
    for (int d=0; d<D; ++d) r[d]=!r[d];
    return r;
  }
  
  vect operator~ () const {
    vect r;
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
    vect<bool,D> r;
    for (int d=0; d<D; ++d) r[d]=(*this)[d]&&x;
    return r;
  }
  
  vect<bool,D> operator|| (const T x) const {
    vect<bool,D> r;
    for (int d=0; d<D; ++d) r[d]=(*this)[d]||x;
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
  
  template<class TT>
  vect<TT,D> ifthen (const vect<TT,D>& a, const vect<TT,D>& b) const {
    vect<TT,D> r;
    for (int d=0; d<D; ++d) r[d]=(*this)[d]?a[d]:b[d];
    return r;
  }
  
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
  
  // Higher order functions
  vect map (T (* const func)(T x)) const {
    vect r;
    for (int d=0; d<D; ++d) r[d] = func(elt[d]);
    return r;
  }
  
  template<class S>
  vect zip (T (* const func)(T x, S y), const vect<S,D>& a) const {
    vect r;
    for (int d=0; d<D; ++d) r[d] = func(elt[d], a[d]);
    return r;
  }
  
  T fold (T (* const func)(T val, T x), T val) const {
    for (int d=0; d<D; ++d) val = func(val, elt[d]);
    return val;
  }
  
  vect scan0 (T (* const func)(T val, T x), T val) const {
    vect r;
    for (int d=0; d<D; ++d) {
      r[d] = val;
      val = func(val, elt[d]);
    }
    return r;
  }
  
  vect scan1 (T (* const func)(T val, T x), T val) const {
    vect r;
    for (int d=0; d<D; ++d) {
      val = func(val, elt[d]);
      r[d] = val;
    }
    return r;
  }
  
  // Input/Output
  void input (istream& is);
  void output (ostream& os) const;
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
  for (int d=1; d<D; ++d) r=max(r,a[d]);
  return r;
}

template<class T,int D>
inline T minval (const vect<T,D>& a) {
  assert (D>0);
  T r(a[0]);
  for (int d=1; d<D; ++d) r=min(r,a[d]);
  return r;
}

template<class T,int D>
inline int maxloc (const vect<T,D>& a) {
  assert (D>0);
  int r(0);
  for (int d=1; d<D; ++d) if (a[d]>a[r]) r=d;
  return r;
}

template<class T,int D>
inline int minloc (const vect<T,D>& a) {
  assert (D>0);
  int r(0);
  for (int d=1; d<D; ++d) if (a[d]<a[r]) r=d;
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

// Higher order functions
template<class T,class TT,int D>
inline vect<TT,D> map (TT (* const func)(T x), const vect<T,D>& a) {
  vect<TT,D> r;
  for (int d=0; d<D; ++d) r[d] = func(a[d]);
  return r;
}
  
template<class T1,class T2,class TT,int D>
inline vect<TT,D> zip (TT (* const func)(T1 x, T2 y),
		       const vect<T1,D>& a, const vect<T2,D>& b)
{
  vect<TT,D> r;
  for (int d=0; d<D; ++d) r[d] = func(a[d], b[d]);
  return r;
}

template<class T,class TT,int D>
inline vect<TT,D> fold (TT (* const func)(TT val, T x), TT val,
			const vect<T,D>& a)
{
  for (int d=0; d<D; ++d) val = func(val, a[d]);
  return val;
}

template<class T,class TT,int D>
inline vect<TT,D> scan0 (TT (* const func)(TT val, T x), TT val,
			 const vect<T,D>& a)
{
  vect<TT,D> r;
  for (int d=0; d<D; ++d) {
    r[d] = val;
    val = func(val, a[d]);
  }
  return r;
}

template<class T,class TT,int D>
inline vect<TT,D> scan1 (TT (* const func)(TT val, T x), TT val,
			 const vect<T,D>& a)
{
  vect<TT,D> r;
  for (int d=0; d<D; ++d) {
    val = func(val, a[d]);
    r[d] = val;
  }
  return r;
}



// Input
template<class T,int D>
inline istream& operator>> (istream& is, vect<T,D>& a) {
  a.input(is);
  return is;
}



// Output
template<class T,int D>
inline ostream& operator<< (ostream& os, const vect<T,D>& a) {
  a.output(os);
  return os;
}



#endif // VECT_HH
