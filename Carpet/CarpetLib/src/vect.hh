#ifndef VECT_HH
#define VECT_HH

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

#include "cctk.h"

using namespace std;



// Forward definition
template<class T, int D> class vect;

// Input/Output
template<class T,int D>
istream& operator>> (istream& is, vect<T,D>& a);
template<class T,int D>
ostream& operator<< (ostream& os, const vect<T,D>& a);



/**
 * A short vector with a size that is specified at compile time.
 */
template<class T, int D>
class vect {
  
  // Fields
  
  /** Vector elements.  */
  T elt[D==0 ? 1 : D];
  
public:
  
  // Constructors
  
  /** Explicit empty constructor.  */
  explicit vect () { }
  
  /** Copy constructor.  */
  vect (const vect& a) {
    for (int d=0; d<D; ++d) elt[d]=a.elt[d];
  }
  
  /** Constructor from a single element.  This constructor might be
      confusing, but it is very convenient.  */
  vect (const T x) {
    for (int d=0; d<D; ++d) elt[d]=x;
  }
  
  /** Constructor for 2-element vectors from 2 elements.  */
  vect (const T x, const T y) {
    assert (D==2);
    // Note: this statement may give "index out of range" warnings.
    // You can safely ignore these.
    elt[0]=x; elt[1]=y;
  }
  
  /** Constructor for 3-element vectors from 3 elements.  */
  vect (const T x, const T y, const T z) {
    assert (D==3);
    // Note: this statement may give "index out of range" warnings.
    // You can safely ignore these.
    elt[0]=x; elt[1]=y; elt[2]=z;
  }
  
  /** Constructor for 4-element vectors from 4 elements.  */
  vect (const T x, const T y, const T z, const T t) {
    assert (D==4);
    // Note: this statement may give "index out of range" warnings.
    // You can safely ignore these.
    elt[0]=x; elt[1]=y; elt[2]=z; elt[3]=t;
  }
  
#if 0
  // This creates confusion
  /** Constructor from a pointer, i.e.\ a C array.  */
  explicit vect (const T* const x) {
    for (int d=0; d<D; ++d) elt[d]=x[d];
  }
#endif
  
  /** Constructor from a vector with a different type.  */
  template<class S>
  /*explicit*/ vect (const vect<S,D>& a) {
    for (int d=0; d<D; ++d) elt[d]=(T)a[d];
  }
  
  /** Create a new 0-element vector with a specific type.  */
  static vect make () {
    assert (D==0);
    return vect();
  }
  
  /** Create a new 1-element vector with a specific type.  */
  static vect make (const T x) {
    assert (D==1);
    return vect(x);
  }
  
  /** Create a new 2-element vector with a specific type.  */
  static vect make (const T x, const T y) {
    assert (D==2);
    return vect(x, y);
  }
  
  /** Create a new 3-element vector with a specific type.  */
  static vect make (const T x, const T y, const T z) {
    assert (D==3);
    return vect(x, y, z);
  }
  
  /** Create a new 4-element vector with a specific type.  */
  static vect make (const T x, const T y, const T z, const T t) {
    assert (D==4);
    return vect(x, y, z, t);
  }
  
  /** Treat a constant pointer as a reference to a constant vector.  */
  static const vect& ref (const T* const x) {
    return *(const vect*)x;
  }
  
  /** Treat a pointer as a reference to a vector.  */
  static vect& ref (T* const x) {
    return *(vect*)x;
  }
  
  /** Create a vector with one element set to 1 and all other elements
      set to zero.  */
  static vect dir (const int d) {
    vect r=(T)0;
    r[d]=1;
    return r;
  }
  
  /** Create a vector with e[i] = i.  */
  static vect seq () {
    vect r;
    for (int d=0; d<D; ++d) r[d]=d;
    return r;
  }
  
  /** Create a vector with e[i] = n + i.  */
  static vect seq (const int n) {
    vect r;
    for (int d=0; d<D; ++d) r[d]=n+d;
    return r;
  }
  
  /** Create a vector with e[i] = n + s * i.  */
  static vect seq (const int n, const int s) {
    vect r;
    for (int d=0; d<D; ++d) r[d]=n+s*d;
    return r;
  }
  
  // Accessors
  
  /** Return a non-writable element of a vector.  */
  // (Don't return a reference; *this might be a temporary)
  // Do return a reference, so that a vector can be accessed as array
  const T& operator[] (const int d) const {
    assert(d>=0 && d<D);
    return elt[d];
  }
  
  /** Return a writable element of a vector as reference.  */
  T& operator[] (const int d) {
    assert(d>=0 && d<D);
    return elt[d];
  }
  
#if 0
  // This creates confusion
  /** Return a pointer to a vector.  */
  operator const T* () const {
    return this;
  }
#endif
  
  /** Return a combination of the vector elements e[a[i]].  The
      element combination is selected by another vector.  */
  template<class TT, int DD>
  vect<T,DD> operator[] (const vect<TT,DD>& a) const {
    vect<T,DD> r;
    // (*this)[] performs index checking
    for (int d=0; d<DD; ++d) r[d] = (*this)[a[d]];
    return r;
  }
  
  // Modifying operators
  vect& operator+=(const T x) {
    for (int d=0; d<D; ++d) elt[d]+=x;
    return *this;
  }
  
  vect& operator-=(const T x) {
    for (int d=0; d<D; ++d) elt[d]-=x;
    return *this;
  }
  
  vect& operator*=(const T x) {
    for (int d=0; d<D; ++d) elt[d]*=x;
    return *this;
  }
  
  vect& operator/=(const T x) {
    for (int d=0; d<D; ++d) elt[d]/=x;
    return *this;
  }
  
  vect& operator%=(const T x) {
    for (int d=0; d<D; ++d) elt[d]%=x;
    return *this;
  }
  
  vect& operator&=(const T x) {
    for (int d=0; d<D; ++d) elt[d]&=x;
    return *this;
  }
  
  vect& operator|=(const T x) {
    for (int d=0; d<D; ++d) elt[d]|=x;
    return *this;
  }
  
  vect& operator^=(const T x) {
    for (int d=0; d<D; ++d) elt[d]^=x;
    return *this;
  }
  
  vect& operator+=(const vect& a) {
    for (int d=0; d<D; ++d) elt[d]+=a[d];
    return *this;
  }
  
  vect& operator-=(const vect& a) {
    for (int d=0; d<D; ++d) elt[d]-=a[d];
    return *this;
  }
  
  vect& operator*=(const vect& a) {
    for (int d=0; d<D; ++d) elt[d]*=a[d];
    return *this;
  }
  
  vect& operator/=(const vect& a) {
    for (int d=0; d<D; ++d) elt[d]/=a[d];
    return *this;
  }
  
  vect& operator%=(const vect& a) {
    for (int d=0; d<D; ++d) elt[d]%=a[d];
    return *this;
  }
  
  vect& operator&=(const vect& a) {
    for (int d=0; d<D; ++d) elt[d]&=a[d];
    return *this;
  }
  
  vect& operator|=(const vect& a) {
    for (int d=0; d<D; ++d) elt[d]|=a[d];
    return *this;
  }
  
  vect& operator^=(const vect& a) {
    for (int d=0; d<D; ++d) elt[d]^=a[d];
    return *this;
  }
  
  // Non-modifying operators
  
  /** Return a new vector where one element has been replaced.  */
  vect replace (const int d, const T x) const {
    assert (d>=0 && d<D);
    vect r;
    for (int dd=0; dd<D; ++dd) r[dd]=dd==d?x:elt[dd];
    return r;
  }
  
  vect operator+ () const {
    vect r;
    for (int d=0; d<D; ++d) r[d]=+elt[d];
    return r;
  }
  
  vect operator- () const {
    vect r;
    for (int d=0; d<D; ++d) r[d]=-elt[d];
    return r;
  }
  
  vect<bool,D> operator! () const {
    vect<bool,D> r;
    for (int d=0; d<D; ++d) r[d]=!elt[d];
    return r;
  }
  
  vect operator~ () const {
    vect r;
    for (int d=0; d<D; ++d) r[d]=~elt[d];
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
    for (int d=0; d<D; ++d) r[d]=elt[d]&&x;
    return r;
  }
  
  vect<bool,D> operator|| (const T x) const {
    vect<bool,D> r;
    for (int d=0; d<D; ++d) r[d]=elt[d]||x;
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
    for (int d=0; d<D; ++d) r[d]=elt[d]&&a[d];
    return r;
  }
  
  vect<bool,D> operator|| (const vect& a) const {
    vect<bool,D> r;
    for (int d=0; d<D; ++d) r[d]=elt[d]||a[d];
    return r;
  }
  
  vect<bool,D> operator== (const vect& a) const {
    vect<bool,D> r;
    for (int d=0; d<D; ++d) r[d]=elt[d]==a[d];
    return r;
  }
  
  vect<bool,D> operator!= (const vect& a) const {
    vect<bool,D> r;
    for (int d=0; d<D; ++d) r[d]=elt[d]!=a[d];
    return r;
  }
  
  vect<bool,D> operator< (const vect& a) const {
    vect<bool,D> r;
    for (int d=0; d<D; ++d) r[d]=elt[d]<a[d];
    return r;
  }
  
  vect<bool,D> operator<= (const vect& a) const {
    vect<bool,D> r;
    for (int d=0; d<D; ++d) r[d]=elt[d]<=a[d];
    return r;
  }
  
  vect<bool,D> operator> (const vect& a) const {
    vect<bool,D> r;
    for (int d=0; d<D; ++d) r[d]=elt[d]>a[d];
    return r;
  }
  
  vect<bool,D> operator>= (const vect& a) const {
    vect<bool,D> r;
    for (int d=0; d<D; ++d) r[d]=elt[d]>=a[d];
    return r;
  }
  
  /** This corresponds to the ?: operator.  Return a vector with the
      elements set to either a[i] or b[i], depending on whether
      (*this)[i] is true or not.  */
  template<class TT>
  vect<TT,D> ifthen (const vect<TT,D>& a, const vect<TT,D>& b) const {
    vect<TT,D> r;
    for (int d=0; d<D; ++d) r[d]=elt[d]?a[d]:b[d];
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
  
  // Input/Output helpers
  void input (istream& is);
  void output (ostream& os) const;
};



// Operators

/** This corresponds to the ?: operator.  Return a vector with the
    elements set to either b[i] or c[i], depending on whether a[i] is
    true or not.  */
template<class T,int D>
inline vect<T,D> either (const vect<bool,D>& a,
                         const vect<T,D>& b, const vect<T,D>& c) {
  vect<T,D> r;
  for (int d=0; d<D; ++d) r[d]=a[d]?b[d]:c[d];
  return r;
}

/** Transpose a vector of a vector */
template<class T, int D, int DD>
inline vect<vect<T,D>,DD> xpose (vect<vect<T,DD>,D> const & a) {
  vect<vect<T,D>,DD> r;
  for (int dd=0; dd<DD; ++dd) for (int d=0; d<D; ++d) r[dd][d] = a[d][dd];
  return r;
}

/** Return the element-wise absolute value.  */
template<class T,int D>
inline vect<T,D> abs (const vect<T,D>& a) {
  vect<T,D> r;
  for (int d=0; d<D; ++d) r[d]=abs(a[d]);
  return r;
}

/** Return the element-wise ceiling.  */
template<class T,int D>
inline vect<T,D> ceil (const vect<T,D>& a) {
  vect<T,D> r;
  for (int d=0; d<D; ++d) r[d]=ceil(a[d]);
  return r;
}

/** Return the element-wise floor.  */
template<class T,int D>
inline vect<T,D> floor (const vect<T,D>& a) {
  vect<T,D> r;
  for (int d=0; d<D; ++d) r[d]=floor(a[d]);
  return r;
}

/** Return the element-wise maximum of two vectors.  */
template<class T,int D>
inline vect<T,D> max (const vect<T,D>& a, const vect<T,D>& b) {
  vect<T,D> r;
  for (int d=0; d<D; ++d) r[d]=max(a[d],b[d]);
  return r;
}

/** Return the element-wise minimum of two vectors.  */
template<class T,int D>
inline vect<T,D> min (const vect<T,D>& a, const vect<T,D>& b) {
  vect<T,D> r;
  for (int d=0; d<D; ++d) r[d]=min(a[d],b[d]);
  return r;
}

/** Return the element-wise power of two vectors.  */
template<class T,class U,int D>
inline vect<T,D> pow (const vect<T,D>& a, const vect<U,D>& b) {
  vect<T,D> r;
  for (int d=0; d<D; ++d) r[d]=pow(a[d],b[d]);
  return r;
}

/** Return the element-wise integer power of two vectors.  */
template<class T,int D>
inline vect<T,D> ipow (const vect<T,D>& a, const vect<int,D>& b) {
  vect<T,D> r;
  for (int d=0; d<D; ++d) r[d]=ipow(a[d],b[d]);
  return r;
}



// Reduction operators

/** Return true iff any of the elements are true (boolean sum).  */
template<int D>
inline bool any (const vect<bool,D>& a) {
  bool r(false);
  for (int d=0; d<D; ++d) r|=a[d];
  return r;
}

/** Return true iff all of the elements are true (boolean product).  */
template<int D>
inline bool all (const vect<bool,D>& a) {
  bool r(true);
  for (int d=0; d<D; ++d) r&=a[d];
  return r;
}

/** Count the number of elements in the vector.  */
template<class T,int D>
inline int count (const vect<T,D>& a) {
  return D;
}

/** Return the dot product of two vectors.  */
template<class T,int D>
inline T dot (const vect<T,D>& a, const vect<T,D>& b) {
  T r(0);
  for (int d=0; d<D; ++d) r+=a[d]*b[d];
  return r;
}

/** Return the Euklidean length.  */
template<class T,int D>
inline T hypot (const vect<T,D>& a) {
  return sqrt(dot(a,a));
}

/** Return the maximum element.  */
template<class T,int D>
inline T maxval (const vect<T,D>& a) {
  assert (D>0);
  T r(a[0]);
  for (int d=1; d<D; ++d) r=max(r,a[d]);
  return r;
}

/** Return the minimum element.  */
template<class T,int D>
inline T minval (const vect<T,D>& a) {
  assert (D>0);
  T r(a[0]);
  for (int d=1; d<D; ++d) r=min(r,a[d]);
  return r;
}

/** Return the index of the first maximum element.  */
template<class T,int D>
inline int maxloc (const vect<T,D>& a) {
  assert (D>0);
  int r(0);
  for (int d=1; d<D; ++d) if (a[d]>a[r]) r=d;
  return r;
}

/** Return the index of the first minimum element.  */
template<class T,int D>
inline int minloc (const vect<T,D>& a) {
  assert (D>0);
  int r(0);
  for (int d=1; d<D; ++d) if (a[d]<a[r]) r=d;
  return r;
}

/** Return the product of the elements.  */
template<class T,int D>
inline T prod (const vect<T,D>& a) {
  T r(1);
  for (int d=0; d<D; ++d) r*=a[d];
  return r;
}

/** Return the size (number of elements) of the vector.  */
template<class T,int D>
inline int size (const vect<T,D>& a) {
  return D;
}

/** Return the sum of the elements.  */
template<class T,int D>
inline T sum (const vect<T,D>& a) {
  T r(0);
  for (int d=0; d<D; ++d) r+=a[d];
  return r;
}

// Higher order functions

/** Return a new vector where the function func() has been applied to
    all elements.  */
template<class T, class U, int D>
inline vect<U,D> map (U (* const func)(T x), const vect<T,D>& a) {
  vect<U,D> r;
  for (int d=0; d<D; ++d) r[d] = func(a[d]);
  return r;
}
  
/** Return a new vector where the function func() has been used
    element-wise to combine a and b.  */
template<class S, class T, class U, int D>
inline vect<U,D> zip (U (* const func)(S x, T y),
                      const vect<S,D>& a, const vect<T,D>& b)
{
  vect<U,D> r;
  for (int d=0; d<D; ++d) r[d] = func(a[d], b[d]);
  return r;
}

/** Return a scalar where the function func() has been used to reduce
    the vector a, starting with the scalar value val.  */
template<class T, class U, int D>
inline U fold (U (* const func)(U val, T x), U val, const vect<T,D>& a)
{
  for (int d=0; d<D; ++d) val = func(val, a[d]);
  return val;
}
  
/** Return a scalar where the function func() has been used to reduce
    the vector a, starting with element 0.  */
template<class T, class U, int D>
inline U fold1 (U (* const func)(U val, T x), const vect<T,D>& a)
{
  assert (D>=1);
  U val = a[0];
  for (int d=1; d<D; ++d) val = func(val, a[d]);
  return val;
}

/** Return a vector where the function func() has been used to scan
    the vector a, starting with the scalar value val.  */
template<class T, class U, int D>
inline vect<U,D> scan0 (U (* const func)(U val, T x), U val,
                        const vect<T,D>& a)
{
  vect<U,D> r;
  for (int d=0; d<D; ++d) {
    r[d] = val;
    val = func(val, a[d]);
  }
  return r;
}

/** Return a vector where the function func() has been used to scan
    the vector a, starting with element 0.  */
template<class T, class U, int D>
inline vect<U,D> scan1 (U (* const func)(U val, T x), U val,
                        const vect<T,D>& a)
{
  vect<U,D> r;
  for (int d=0; d<D; ++d) {
    val = func(val, a[d]);
    r[d] = val;
  }
  return r;
}



// Input

/** Read a formatted vector from a stream.  */
template<class T,int D>
inline istream& operator>> (istream& is, vect<T,D>& a) {
  a.input(is);
  return is;
}



// Output

/** Write a vector formatted to a stream.  */
template<class T,int D>
inline ostream& operator<< (ostream& os, const vect<T,D>& a) {
  a.output(os);
  return os;
}



#if 0  
// Specialise explicit constructors

/** Constructor for 2-element vectors from 2 elements.  */
template<class T>
inline vect<T,2>::vect<T,2> (const T x, const T y) {
  elt[0]=x; elt[1]=y;
}

/** Constructor for 3-element vectors from 3 elements.  */
vect (const T x, const T y, const T z) {
  assert (D==3);
  elt[0]=x; elt[1]=y; elt[2]=z;
}

/** Constructor for 4-element vectors from 4 elements.  */
vect (const T x, const T y, const T z, const T t) {
  assert (D==4);
  elt[0]=x; elt[1]=y; elt[2]=z; elt[3]=t;
}
#endif



// Specialise some constructors for lower dimensions
// These functions are declared, but never defined, so that using them
// will result in a linker error

template<> vect<int,0>::vect (const int x, const int y);
template<> vect<int,1>::vect (const int x, const int y);

template<> vect<int,0>::vect (const int x, const int y, const int z);
template<> vect<int,1>::vect (const int x, const int y, const int z);
template<> vect<int,2>::vect (const int x, const int y, const int z);

template<> vect<int,0>::vect (const int x, const int y, const int z, const int t);
template<> vect<int,1>::vect (const int x, const int y, const int z, const int t);
template<> vect<int,2>::vect (const int x, const int y, const int z, const int t);
template<> vect<int,3>::vect (const int x, const int y, const int z, const int t);




// Specialise for CCTK_REAL

template<>
inline vect<CCTK_REAL,3>& vect<CCTK_REAL,3>::operator%=(const vect<CCTK_REAL,3>& a) {
  for (int d=0; d<3; ++d) {
    elt[d]=fmod(elt[d],a[d]);
    if (elt[d]>a[d]*(CCTK_REAL)(1.0-1.0e-10)) elt[d]=(CCTK_REAL)0;
    if (elt[d]<a[d]*(CCTK_REAL)(    1.0e-10)) elt[d]=(CCTK_REAL)0;
  }
  return *this;
}



#endif // VECT_HH
