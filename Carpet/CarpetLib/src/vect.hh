#ifndef VECT_HH
#define VECT_HH

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

#include "cctk.h"

#include "vect_helpers.hh"

using namespace std;



#if 0

// A pure function returns a value that depends only on the function
// arguments and on global variables, and the function has no side
// effects.
#ifdef HAVE_CCTK_CXX_ATTRIBUTE_PURE
#  define PURE __attribute__((pure))
#else
#  define PURE
#endif

// A const function returns a value that depends only on the function
// arguments, and the function has no side effects.  This is even more
// strict than pure functions.  Const functions cannot dereference
// pointers or references (or this).
#ifdef HAVE_CCTK_CXX_ATTRIBUTE_CONST
#  define CONST __attribute__((const))
#else
#  define CONST
#endif

#else

// Don't take any risks
#  define PURE
#  define CONST

#endif



// Forward definition
template<typename T, int D> class vect;

// Input/Output
template<typename T,int D>
istream& operator>> (istream& is, vect<T,D>& a);
template<typename T,int D>
ostream& operator<< (ostream& os, const vect<T,D>& a);



/**
 * A short vector with a size that is specified at compile time.
 */
template<typename T, int D>
class vect {
  
  // Fields
  
  /** Vector elements.  */
  T elt[D==0 ? 1 : D];
  
public:
  
  // Constructors
  
  /** Explicit empty constructor.  */
  explicit vect () CONST { }
  
  /** Copy constructor.  */
  vect (const vect& a) PURE {
    for (int d=0; d<D; ++d) elt[d]=a.elt[d];
  }
  
  /** Constructor from a single element.  This constructor might be
      confusing, but it is very convenient.  */
  vect (const T& x) PURE {
    for (int d=0; d<D; ++d) elt[d]=x;
  }
  
  /** Constructor for 2-element vectors from 2 elements.  */
  vect (const T& x, const T& y) PURE {
    assert (D==2);
    // Note: this statement may give "index out of range" warnings.
    // You can safely ignore these.
    elt[0]=x; elt[1]=y;
  }
  
  /** Constructor for 3-element vectors from 3 elements.  */
  vect (const T& x, const T& y, const T& z) PURE {
    assert (D==3);
    // Note: this statement may give "index out of range" warnings.
    // You can safely ignore these.
    elt[0]=x; elt[1]=y; elt[2]=z;
  }
  
  /** Constructor for 4-element vectors from 4 elements.  */
  vect (const T& x, const T& y, const T& z, const T& t) PURE {
    assert (D==4);
    // Note: this statement may give "index out of range" warnings.
    // You can safely ignore these.
    elt[0]=x; elt[1]=y; elt[2]=z; elt[3]=t;
  }
  
#if 0
  // This creates confusion
  /** Constructor from a pointer, i.e., a C array.  */
  explicit vect (const T* const x) PURE {
    for (int d=0; d<D; ++d) elt[d]=x[d];
  }
#endif
  
#if 0
  // This leads to an ICE on AIX
  template<int E>
  operator vect<vect<T,D>,E> () CONST {
    vect<vect<T,D>,E> r;
    for (int e=0; e<E; ++e) r[e]=*this;
    return r;
  }
#endif
  
  /** Constructor from a vector with a different type.  */
  template<typename S>
  /*explicit*/ vect (const vect<S,D>& a) /*PURE*/ {
    for (int d=0; d<D; ++d) elt[d]=(T)a[d];
  }
  
  /** Create a new 0-element vector with a specific type.  */
  static vect make () CONST {
    assert (D==0);
    return vect();
  }
  
  /** Create a new 1-element vector with a specific type.  */
  static vect make (const T& x) PURE {
    assert (D==1);
    return vect(x);
  }
  
  /** Create a new 2-element vector with a specific type.  */
  static vect make (const T& x, const T& y) PURE {
    assert (D==2);
    return vect(x, y);
  }
  
  /** Create a new 3-element vector with a specific type.  */
  static vect make (const T& x, const T& y, const T& z) PURE {
    assert (D==3);
    return vect(x, y, z);
  }
  
  /** Create a new 4-element vector with a specific type.  */
  static vect make (const T& x, const T& y, const T& z, const T& t) PURE {
    assert (D==4);
    return vect(x, y, z, t);
  }
  
  /** Treat a constant pointer as a reference to a constant vector.  */
  static const vect& ref (const T* const x) PURE {
    return *(const vect*)x;
  }
  
  /** Treat a pointer as a reference to a vector.  */
  static vect& ref (T* const x) PURE {
    return *(vect*)x;
  }
  
  /** Create a vector with one element set to 1 and all other elements
      set to zero.  */
  static vect dir (const int d) CONST {
    vect r=(T)0;
    r[d]=1;
    return r;
  }
  
  /** Create a vector with e[i] = i.  */
  static vect seq () CONST {
    vect r;
    for (int d=0; d<D; ++d) r[d]=d;
    return r;
  }
  
  /** Create a vector with e[i] = n + i.  */
  static vect seq (const int n) CONST {
    vect r;
    for (int d=0; d<D; ++d) r[d]=n+d;
    return r;
  }
  
  /** Create a vector with e[i] = n + s * i.  */
  static vect seq (const int n, const int s) CONST {
    vect r;
    for (int d=0; d<D; ++d) r[d]=n+s*d;
    return r;
  }
  
  // Accessors
  
  /** Return a non-writable element of a vector.  */
  // (Don't return a reference; *this might be a temporary)
  // Do return a reference, so that a vector can be accessed as array
  const T& operator[] (const int d) const PURE {
    assert(d>=0 && d<D);
    return elt[d];
  }
  
  /** Return a writable element of a vector as reference.  */
  T& operator[] (const int d) PURE {
    assert(d>=0 && d<D);
    return elt[d];
  }
  
#if 0
  // This creates confusion
  /** Return a pointer to a vector.  */
  operator const T* () const PURE {
    return this;
  }
#endif
  
  /** Return a combination of the vector elements e[a[i]].  The
      element combination is selected by another vector.  */
  template<typename TT, int DD>
  vect<T,DD> operator[] (const vect<TT,DD>& a) const /*PURE*/ {
    vect<T,DD> r;
    // (*this)[] performs index checking
    for (int d=0; d<DD; ++d) r[d] = (*this)[a[d]];
    return r;
  }
  
  // Modifying operators
  DECLARE_MEMBER_OPERATOR_1_REF (operator+=, +=);
  DECLARE_MEMBER_OPERATOR_1_REF (operator-=, -=);
  DECLARE_MEMBER_OPERATOR_1_REF (operator*=, *=);
  DECLARE_MEMBER_OPERATOR_1_REF (operator/=, /=);
  DECLARE_MEMBER_OPERATOR_1_REF (operator%=, %=);
  DECLARE_MEMBER_OPERATOR_1_REF (operator&=, &=);
  DECLARE_MEMBER_OPERATOR_1_REF (operator|=, |=);
  DECLARE_MEMBER_OPERATOR_1_REF (operator^=, ^=);
  
  // Non-modifying operators
  
  /** Return a new vector where one element has been replaced.  */
  vect replace (const int d, const T& x) const PURE {
    assert (d>=0 && d<D);
    vect r;
    for (int dd=0; dd<D; ++dd) r[dd]=dd==d?x:elt[dd];
    return r;
  }
  
  vect reverse () const PURE {
    vect r;
    for (int d=0; d<D; ++d) r[d]=elt[D-1-d];
    return r;
  }
  
  DECLARE_MEMBER_OPERATOR_0 (operator+, +)
  DECLARE_MEMBER_OPERATOR_0 (operator-, -)
  DECLARE_MEMBER_OPERATOR_0 (operator~, ~)
  // DECLARE_MEMBER_OPERATOR_0_RET (operator!, !, bool)
  
#if 0
  /** This corresponds to the ?: operator.  Return a vector with the
      elements set to either a[i] or b[i], depending on whether
      (*this)[i] is true or not.  */
  template<typename TT>
  vect<TT,D> ifthen (const vect<TT,D>& a, const vect<TT,D>& b) const /*PURE*/ {
    vect<TT,D> r;
    for (int d=0; d<D; ++d) r[d]=elt[d]?a[d]:b[d];
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
    iter (vect &a) PURE: vec(a), d(0) { }
    iter& operator++ () { assert(d<D); ++d; return *this; }
    bool operator bool () const PURE { return d==D; }
    T& operator* () PURE { return vec[d]; }
  };
#endif
  
  // Memory usage
  size_t memory () const { return D * memoryof (*elt); }
  
  // Input/Output helpers
  void input (istream& is);
  void output (ostream& os) const;
};



// Operators

/** This corresponds to the ?: operator.  Return a vector with the
    elements set to either b[i] or c[i], depending on whether a[i] is
    true or not.  */
template<typename S,typename T,int D>
inline vect<T,D> either (const vect<S,D>& a,
                         const vect<T,D>& b, const vect<T,D>& c) PURE;
template<typename S,typename T,int D>
inline vect<T,D> either (const vect<S,D>& a,
                         const vect<T,D>& b, const vect<T,D>& c)
{
  vect<T,D> r;
  for (int d=0; d<D; ++d) r[d]=a[d]?b[d]:c[d];
  return r;
}

template<typename S,typename T,int D>
inline vect<T,D> either (const vect<S,D>& a,
                         const T& b, const T& c) PURE;
template<typename S,typename T,int D>
inline vect<T,D> either (const vect<S,D>& a,
                         const T& b, const T& c)
{
  return either (a, vect<T,D>(b), vect<T,D>(c));
}

/** Transpose a vector of a vector */
template<typename T, int D, int DD>
inline vect<vect<T,D>,DD> xpose (vect<vect<T,DD>,D> const & a) PURE;
template<typename T, int D, int DD>
inline vect<vect<T,D>,DD> xpose (vect<vect<T,DD>,D> const & a) {
  vect<vect<T,D>,DD> r;
  for (int dd=0; dd<DD; ++dd) for (int d=0; d<D; ++d) r[dd][d] = a[d][dd];
  return r;
}

/** Return the element-wise integer power of two vectors.  */
template<typename T,int D>
inline vect<T,D> ipow (const vect<T,D>& a, const vect<int,D>& b) PURE;
template<typename T,int D>
inline vect<T,D> ipow (const vect<T,D>& a, const vect<int,D>& b) {
  vect<T,D> r;
  for (int d=0; d<D; ++d) r[d]=ipow(a[d],b[d]);
  return r;
}

DECLARE_FUNCTION_1 (abs)
DECLARE_FUNCTION_1 (ceil)
DECLARE_FUNCTION_1 (floor)

DECLARE_OPERATOR_1_RET (operator!, !, bool)

DECLARE_FUNCTION_2 (max)
DECLARE_FUNCTION_2 (min)
DECLARE_FUNCTION_2 (pow)

DECLARE_OPERATOR_2 (operator+, +)
DECLARE_OPERATOR_2 (operator-, -)
DECLARE_OPERATOR_2 (operator*, *)
DECLARE_OPERATOR_2 (operator/, /)
DECLARE_OPERATOR_2 (operator%, %)
DECLARE_OPERATOR_2 (operator&, &)
DECLARE_OPERATOR_2 (operator|, |)
DECLARE_OPERATOR_2 (operator^, ^)

DECLARE_OPERATOR_2_RET (operator&&, &&, bool)
DECLARE_OPERATOR_2_RET (operator||, ||, bool)
DECLARE_OPERATOR_2_RET (operator==, ==, bool)
DECLARE_OPERATOR_2_RET (operator!=, !=, bool)
DECLARE_OPERATOR_2_RET (operator< , < , bool)
DECLARE_OPERATOR_2_RET (operator<=, <=, bool)
DECLARE_OPERATOR_2_RET (operator> , > , bool)
DECLARE_OPERATOR_2_RET (operator>=, >=, bool)



// Reduction operators

// Identity
#define id(x) (x)

DECLARE_REDUCTION_OPERATOR_1_T_RET (all,true ,&=,id,bool,bool)
DECLARE_REDUCTION_OPERATOR_1_T_RET (any,false,|=,id,bool,bool)

DECLARE_REDUCTION_FUNCTION_1 (maxval,a[0],max,id)
DECLARE_REDUCTION_FUNCTION_1 (minval,a[0],min,id)
DECLARE_REDUCTION_OPERATOR_1 (prod,1,*=,id)
DECLARE_REDUCTION_OPERATOR_1 (sum,0,+=,id)

DECLARE_REDUCTION_OPERATOR_2 (dot  ,0,+=,*,id  )
DECLARE_REDUCTION_OPERATOR_2 (hypot,0,+=,*,sqrt)

/** Count the number of elements in the vector.  */
template<typename T,int D>
inline int count (const vect<T,D>& a) PURE;
template<typename T,int D>
inline int count (const vect<T,D>& a) {
  return D;
}

/** Return the size (number of elements) of the vector.  */
template<typename T,int D>
inline int size (const vect<T,D>& a) CONST;
template<typename T,int D>
inline int size (const vect<T,D>& a) {
  return D;
}

/** Return the index of the first maximum element.  */
template<typename T,int D>
inline int maxloc (const vect<T,D>& a) PURE;
template<typename T,int D>
inline int maxloc (const vect<T,D>& a) {
  assert (D>0);
  int r(0);
  for (int d=1; d<D; ++d) if (a[d]>a[r]) r=d;
  return r;
}

/** Return the index of the first minimum element.  */
template<typename T,int D>
inline int minloc (const vect<T,D>& a) PURE;
template<typename T,int D>
inline int minloc (const vect<T,D>& a) {
  assert (D>0);
  int r(0);
  for (int d=1; d<D; ++d) if (a[d]<a[r]) r=d;
  return r;
}

/** Return the n-dimensional linear array index.  */
template<typename T,int D>
inline T index (const vect<T,D>& lsh, const vect<T,D>& ind) PURE;
template<typename T,int D>
inline T index (const vect<T,D>& lsh, const vect<T,D>& ind) {
  T r(0);
  for (int d=D-1; d>=0; --d) {
    assert (lsh[d]>=0);
    // Be generous, and allow relative indices which may be negtive
    // assert (ind[d]>=0 and ind[d]<lsh[d]);
    assert (abs(ind[d])<=lsh[d]);
    r = r * lsh[d] + ind[d];
  }
  return r;
}



// Higher order functions

#if 0
// They are rarely used, so disable them

/** Return a new vector where the function func() has been applied to
    all elements.  */
template<typename T, typename U, int D>
inline vect<U,D> map (U (* const func)(T x), const vect<T,D>& a) {
  vect<U,D> r;
  for (int d=0; d<D; ++d) r[d] = func(a[d]);
  return r;
}
  
/** Return a new vector where the function func() has been used
    element-wise to combine a and b.  */
template<typename S, typename T, typename U, int D>
inline vect<U,D> zip (U (* const func)(S x, T y),
                      const vect<S,D>& a, const vect<T,D>& b)
{
  vect<U,D> r;
  for (int d=0; d<D; ++d) r[d] = func(a[d], b[d]);
  return r;
}

/** Return a scalar where the function func() has been used to reduce
    the vector a, starting with the scalar value val.  */
template<typename T, typename U, int D>
inline U fold (U (* const func)(U val, T x), U val, const vect<T,D>& a)
{
  for (int d=0; d<D; ++d) val = func(val, a[d]);
  return val;
}
  
/** Return a scalar where the function func() has been used to reduce
    the vector a, starting with element 0.  */
template<typename T, typename U, int D>
inline U fold1 (U (* const func)(U val, T x), const vect<T,D>& a)
{
  assert (D>=1);
  U val = a[0];
  for (int d=1; d<D; ++d) val = func(val, a[d]);
  return val;
}

/** Return a vector where the function func() has been used to scan
    the vector a, starting with the scalar value val.  */
template<typename T, typename U, int D>
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
template<typename T, typename U, int D>
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
#endif



// Memory usage

template<typename T,int D>
inline size_t memoryof (vect<T,D> const & a) { return a.memory(); }



// Input

/** Read a formatted vector from a stream.  */
template<typename T,int D>
inline istream& operator>> (istream& is, vect<T,D>& a) {
  a.input(is);
  return is;
}



// Output

/** Write a vector formatted to a stream.  */
template<typename T,int D>
inline ostream& operator<< (ostream& os, const vect<T,D>& a) {
  a.output(os);
  return os;
}



#if 0  
// Specialise explicit constructors

/** Constructor for 2-element vectors from 2 elements.  */
template<typename T>
inline vect<T,2>::vect<T,2> (const T& x, const T& y) PURE;
template<typename T>
inline vect<T,2>::vect<T,2> (const T& x, const T& y) {
  elt[0]=x; elt[1]=y;
}

/** Constructor for 3-element vectors from 3 elements.  */
vect (const T& x, const T& y, const T& z) PURE;
vect (const T& x, const T& y, const T& z) {
  assert (D==3);
  elt[0]=x; elt[1]=y; elt[2]=z;
}

/** Constructor for 4-element vectors from 4 elements.  */
vect (const T& x, const T& y, const T& z, const T& t) PURE;
vect (const T& x, const T& y, const T& z, const T& t) {
  assert (D==4);
  elt[0]=x; elt[1]=y; elt[2]=z; elt[3]=t;
}
#endif



////////////////////////////////////////////////////////////////////////////////



// Specialise some constructors for lower dimensions
// These functions are declared, but never defined, so that using them
// will result in a linker error

template<> inline vect<int,0>::vect (const int& x, const int& y) { assert(0); }
template<> inline vect<int,1>::vect (const int& x, const int& y) { assert(0); }

template<> inline vect<int,0>::vect (const int& x, const int& y, const int& z) { assert(0); }
template<> inline vect<int,1>::vect (const int& x, const int& y, const int& z) { assert(0); }
template<> inline vect<int,2>::vect (const int& x, const int& y, const int& z) { assert(0); }

template<> inline vect<int,0>::vect (const int& x, const int& y, const int& z, const int& t) { assert(0); }
template<> inline vect<int,1>::vect (const int& x, const int& y, const int& z, const int& t) { assert(0); }
template<> inline vect<int,2>::vect (const int& x, const int& y, const int& z, const int& t) { assert(0); }
template<> inline vect<int,3>::vect (const int& x, const int& y, const int& z, const int& t) { assert(0); }




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

template<>
inline vect<CCTK_REAL,3> operator%(const vect<CCTK_REAL,3>& a, const vect<CCTK_REAL,3>& b) {
  vect<CCTK_REAL,3> r;
  for (int d=0; d<3; ++d) {
    r[d]=fmod(a[d],b[d]);
    if (r[d]>b[d]*(CCTK_REAL)(1.0-1.0e-10)) r[d]=(CCTK_REAL)0;
    if (r[d]<b[d]*(CCTK_REAL)(    1.0e-10)) r[d]=(CCTK_REAL)0;
  }
  return r;
}



#endif // VECT_HH
