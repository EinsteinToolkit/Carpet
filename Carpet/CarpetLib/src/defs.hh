#ifndef DEFS_HH
#define DEFS_HH

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <stack>
#include <vector>

#include "cctk.h"

#include "typeprops.hh"



using namespace std;



// Stringify
#define STRINGIFY1(x) #x
#define STRINGIFY(x) STRINGIFY1(x)



// Define the restrict qualifier
#ifdef CCTK_CXX_RESTRICT
#  define restrict CCTK_CXX_RESTRICT
#endif



// Structure member offsets
#undef offsetof
#define offsetof(TYPE,MEMBER) ((size_t)&((TYPE*)0)->MEMBER)
#undef __offsetof__
#define __offsetof__ offsetof



// Number of dimensions
#ifndef CARPET_DIM
#  define CARPET_DIM 3
#endif
const int dim = CARPET_DIM;



// Begin a new line without flushing the output buffer
char const * const eol = "\n";


  
// A compile time pseudo assert statement
#define static_assert(_x, _msg) do { typedef int ai[(_x) ? 1 : -1]; } while(0)



// Check a return value
#define check(_expr) do { bool const _val = (_expr); assert(_val); } while(0)



// Use this macro AT instead of vector's operator[] or at().
// Depending on the macro CARPET_OPTIMISE, this macro AT either checks
// for valid indices or not.
#if ! defined(CARPET_OPTIMISE)
#  define AT(index) at(index)
#else
#  define AT(index) operator[](index)
#endif
  


// Some shortcuts for type names
template<typename T, int D> class vect;
template<typename T, int D> class bbox;
template<typename T, int D> class bboxset;
template<typename T, int D, typename P> class fulltree;

typedef vect<bool,dim>      bvect;
typedef vect<int,dim>       ivect;
typedef vect<CCTK_INT,dim>  jvect;
typedef vect<CCTK_REAL,dim> rvect;
typedef bbox<int,dim>       ibbox;
typedef bbox<CCTK_INT,dim>  jbbox;
typedef bbox<CCTK_REAL,dim> rbbox;
typedef bboxset<int,dim>    ibset;
  
// (Try to replace these by b2vect and i2vect)
typedef vect<vect<bool,2>,dim> bbvect;
typedef vect<vect<int,2>,dim>  iivect;
typedef vect<vect<CCTK_INT,2>,dim> jjvect;

typedef vect<vect<bool,dim>,2> b2vect;
typedef vect<vect<int,dim>,2>  i2vect;



struct pseudoregion_t;
struct region_t;

typedef fulltree<int,dim,pseudoregion_t> ipfulltree;



// A general type
enum centering { error_centered, vertex_centered, cell_centered };



// Useful helper
template<class T>
inline T square (const T& x) CCTK_ATTRIBUTE_CONST;
template<class T>
inline T square (const T& x) { return x*x; }

// Another useful helper
template<class T>
T ipow (T x, int y) CCTK_ATTRIBUTE_CONST;



// Access to CarpetLib parameters
CCTK_INT get_poison_value() CCTK_ATTRIBUTE_CONST;
CCTK_INT get_deadbeef() CCTK_ATTRIBUTE_CONST;



// Input streams
struct input_error { };
void skipws (istream& is);
void expect (istream& is, char c);
void consume (istream& is, char c);
void consume (istream& is, char const * c);



// Names for types

#ifdef HAVE_CCTK_INT1
inline const char * typestring (const CCTK_INT1&) CCTK_ATTRIBUTE_CONST;
inline const char * typestring (const CCTK_INT1&)
{ return "CCTK_INT1"; }
#endif

#ifdef HAVE_CCTK_INT2
inline const char * typestring (const CCTK_INT2&) CCTK_ATTRIBUTE_CONST;
inline const char * typestring (const CCTK_INT2&)
{ return "CCTK_INT2"; }
#endif

#ifdef HAVE_CCTK_INT4
inline const char * typestring (const CCTK_INT4&) CCTK_ATTRIBUTE_CONST;
inline const char * typestring (const CCTK_INT4&)
{ return "CCTK_INT4"; }
#endif

#ifdef HAVE_CCTK_INT8
inline const char * typestring (const CCTK_INT8&) CCTK_ATTRIBUTE_CONST;
inline const char * typestring (const CCTK_INT8&)
{ return "CCTK_INT8"; }
#endif

#ifdef HAVE_CCTK_REAL4
inline const char * typestring (const CCTK_REAL4&) CCTK_ATTRIBUTE_CONST;
inline const char * typestring (const CCTK_REAL4&)
{ return "CCTK_REAL4"; }
#endif

#ifdef HAVE_CCTK_REAL8
inline const char * typestring (const CCTK_REAL8&) CCTK_ATTRIBUTE_CONST;
inline const char * typestring (const CCTK_REAL8&)
{ return "CCTK_REAL8"; }
#endif

#ifdef HAVE_CCTK_REAL16
inline const char * typestring (const CCTK_REAL16&) CCTK_ATTRIBUTE_CONST;
inline const char * typestring (const CCTK_REAL16&)
{ return "CCTK_REAL16"; }
#endif

#ifdef HAVE_CCTK_REAL4
inline const char * typestring (const CCTK_COMPLEX8&) CCTK_ATTRIBUTE_CONST;
inline const char * typestring (const CCTK_COMPLEX8&)
{ return "CCTK_COMPLEX8"; }
#endif

#ifdef HAVE_CCTK_REAL8
inline const char * typestring (const CCTK_COMPLEX16&) CCTK_ATTRIBUTE_CONST;
inline const char * typestring (const CCTK_COMPLEX16&)
{ return "CCTK_COMPLEX16"; }
#endif

#ifdef HAVE_CCTK_REAL16
inline const char * typestring (const CCTK_COMPLEX32&) CCTK_ATTRIBUTE_CONST;
inline const char * typestring (const CCTK_COMPLEX32&)
{ return "CCTK_COMPLEX32"; }
#endif



// Capture the system's isnan function

#ifdef HAVE_CCTK_REAL4
inline int myisnan (CCTK_REAL4 const & x) CCTK_ATTRIBUTE_CONST;
int myisnan (CCTK_REAL4 const & x)
{ return isnan (x); }
#endif
#ifdef HAVE_CCTK_REAL8
inline int myisnan (CCTK_REAL8 const & x) CCTK_ATTRIBUTE_CONST;
inline int myisnan (CCTK_REAL8 const & x)
{ return isnan (x); }
#endif
#ifdef HAVE_CCTK_REAL16
inline int myisnan (CCTK_REAL16 const & x) CCTK_ATTRIBUTE_CONST;
inline int myisnan (CCTK_REAL16 const & x)
{ return isnan (x); }
#endif

#undef isnan



namespace CarpetLib {
  namespace good {
    
    // Explicitly overload some functions for all types in the same
    // namespace CarpetLib::good, to circumvent confusion among some
    // compilers
    
    //
    // abs
    //
    
    template <typename T>
    inline typename typeprops<T>::real abs (T const & x) CCTK_ATTRIBUTE_CONST;
    template <typename T>
    inline typename typeprops<T>::real abs (T const & x)
    { return std::abs (x); }
    
//     // This does not work on Linux with Intel compilers, which do not
//     // always have long long llabs (long long)
//     template<> inline signed char abs<signed char> (signed char const & x) CCTK_ATTRIBUTE_CONST { return ::abs (x); }
//     template<> inline unsigned char abs<unsigned char> (unsigned char const & x) CCTK_ATTRIBUTE_CONST { return ::abs (x); }
//     template<> inline short abs<short> (short const & x) { return ::abs (x); }
//     template<> inline int abs<int> (int const & x) CCTK_ATTRIBUTE_CONST { return ::abs (x); }
//     template<> inline long abs<long> (long const & x) CCTK_ATTRIBUTE_CONST { return ::labs (x); }
// #ifdef SIZEOF_LONG_LONG
//     inline long long abs<long long> (long long const & x) CCTK_ATTRIBUTE_CONST { return ::llabs (x); }
// #endif
    
//     // This template does not work on AIX, which does not have long
//     // long abs (long long)
// #ifdef HAVE_CCTK_INT1
//     template<> inline CCTK_INT1 abs<CCTK_INT1> (CCTK_INT1 const & x) CCTK_ATTRIBUTE_CONST { return x < 0 ? - x : x; }
// #endif
// #ifdef HAVE_CCTK_INT2
//     template<> inline CCTK_INT2 abs<CCTK_INT2> (CCTK_INT2 const & x) CCTK_ATTRIBUTE_CONST { return x < 0 ? - x : x; }
// #endif
// #ifdef HAVE_CCTK_INT4
//     template<> inline CCTK_INT4 abs<CCTK_INT4> (CCTK_INT4 const & x) CCTK_ATTRIBUTE_CONST { return x < 0 ? - x : x; }
// #endif
// #ifdef HAVE_CCTK_INT8
//     template<> inline CCTK_INT8 abs<CCTK_INT8> (CCTK_INT8 const & x) CCTK_ATTRIBUTE_CONST { return x < 0 ? - x : x; }
// #endif
    
#ifdef HAVE_CCTK_COMPLEX8
    template<> inline CCTK_REAL4 abs<CCTK_COMPLEX8> (CCTK_COMPLEX8 const & x) CCTK_ATTRIBUTE_CONST;
    template<> inline CCTK_REAL4 abs<CCTK_COMPLEX8> (CCTK_COMPLEX8 const & x)
    { return CCTK_Cmplx8Abs (x); }
#endif
#ifdef HAVE_CCTK_COMPLEX16
    template<> inline CCTK_REAL8 abs<CCTK_COMPLEX16> (CCTK_COMPLEX16 const & x) CCTK_ATTRIBUTE_CONST;
    template<> inline CCTK_REAL8 abs<CCTK_COMPLEX16> (CCTK_COMPLEX16 const & x)
    { return CCTK_Cmplx16Abs (x); }
#endif
#ifdef HAVE_CCTK_COMPLEX32
    template<> inline CCTK_REAL16 abs<CCTK_COMPLEX32> (CCTK_COMPLEX32 const & x) CCTK_ATTRIBUTE_CONST;
    template<> inline CCTK_REAL16 abs<CCTK_COMPLEX32> (CCTK_COMPLEX32 const & x)
    { return CCTK_Cmplx32Abs (x); }
#endif
    
    //
    // isnan
    //
    
    // Default implementation, only good for integers
    template <typename T>
    inline int isnan (T const & x) CCTK_ATTRIBUTE_CONST;
    template <typename T>
    inline int isnan (T const & x)
    { return 0; }
    
#ifdef HAVE_CCTK_REAL4
    template<> inline int isnan (CCTK_REAL4 const & x) CCTK_ATTRIBUTE_CONST;
    template<> inline int isnan (CCTK_REAL4 const & x)
    { return myisnan (x); }
#endif
#ifdef HAVE_CCTK_REAL8
    template<> inline int isnan (CCTK_REAL8 const & x) CCTK_ATTRIBUTE_CONST;
    template<> inline int isnan (CCTK_REAL8 const & x)
    { return myisnan (x); }
#endif
#ifdef HAVE_CCTK_REAL16
    template<> inline int isnan (CCTK_REAL16 const & x) CCTK_ATTRIBUTE_CONST;
    template<> inline int isnan (CCTK_REAL16 const & x)
    { return myisnan (x); }
#endif
    
#ifdef HAVE_CCTK_COMPLEX8
    template<> inline int isnan (CCTK_COMPLEX8 const & x) CCTK_ATTRIBUTE_CONST;
    template<> inline int isnan (CCTK_COMPLEX8 const & x)
    { return myisnan (CCTK_Cmplx8Real (x)) or myisnan (CCTK_Cmplx8Imag (x)); }
#endif
#ifdef HAVE_CCTK_COMPLEX16
    template<> inline int isnan (CCTK_COMPLEX16 const & x) CCTK_ATTRIBUTE_CONST;
    template<> inline int isnan (CCTK_COMPLEX16 const & x)
    { return myisnan (CCTK_Cmplx16Real (x)) or myisnan (CCTK_Cmplx16Imag (x)); }
#endif
#ifdef HAVE_CCTK_COMPLEX32
    template<> inline int isnan (CCTK_COMPLEX32 const & x) CCTK_ATTRIBUTE_CONST;
    template<> inline int isnan (CCTK_COMPLEX32 const & x)
    { return myisnan (CCTK_Cmplx32Real (x)) or myisnan (CCTK_Cmplx32Imag (x)); }
#endif
    
  } // namespace good
} // namespace CarpetLib



// Container memory usage
inline size_t memoryof (char const e) CCTK_ATTRIBUTE_CONST;
inline size_t memoryof (short const e) CCTK_ATTRIBUTE_CONST;
inline size_t memoryof (int const e) CCTK_ATTRIBUTE_CONST;
inline size_t memoryof (long const e) CCTK_ATTRIBUTE_CONST;
inline size_t memoryof (long long const e) CCTK_ATTRIBUTE_CONST;
inline size_t memoryof (unsigned char const e) CCTK_ATTRIBUTE_CONST;
inline size_t memoryof (unsigned short const e) CCTK_ATTRIBUTE_CONST;
inline size_t memoryof (unsigned int const e) CCTK_ATTRIBUTE_CONST;
inline size_t memoryof (unsigned long const e) CCTK_ATTRIBUTE_CONST;
inline size_t memoryof (unsigned long long const e) CCTK_ATTRIBUTE_CONST;
inline size_t memoryof (float const e) CCTK_ATTRIBUTE_CONST;
inline size_t memoryof (double const e) CCTK_ATTRIBUTE_CONST;
inline size_t memoryof (long double const e) CCTK_ATTRIBUTE_CONST;
inline size_t memoryof (void * const e) CCTK_ATTRIBUTE_CONST;
inline size_t memoryof (void const * const e) CCTK_ATTRIBUTE_CONST;
template<class T> inline size_t memoryof (T * const e) CCTK_ATTRIBUTE_CONST;
template<class T> inline size_t memoryof (T const * const e) CCTK_ATTRIBUTE_CONST;
template<class T> inline size_t memoryof (typename list<T>::iterator const & i) CCTK_ATTRIBUTE_CONST;
template<class T> inline size_t memoryof (typename list<T>::const_iterator const & i) CCTK_ATTRIBUTE_CONST;

inline size_t memoryof (char const e) { return sizeof e; }
inline size_t memoryof (short const e) { return sizeof e; }
inline size_t memoryof (int const e) { return sizeof e; }
inline size_t memoryof (long const e) { return sizeof e; }
inline size_t memoryof (long long const e) { return sizeof e; }
inline size_t memoryof (unsigned char const e) { return sizeof e; }
inline size_t memoryof (unsigned short const e) { return sizeof e; }
inline size_t memoryof (unsigned int const e) { return sizeof e; }
inline size_t memoryof (unsigned long const e) { return sizeof e; }
inline size_t memoryof (unsigned long long const e) { return sizeof e; }
inline size_t memoryof (float const e) { return sizeof e; }
inline size_t memoryof (double const e) { return sizeof e; }
inline size_t memoryof (long double const e) { return sizeof e; }
inline size_t memoryof (void * const e) { return sizeof e; }
inline size_t memoryof (void const * const e) { return sizeof e; }
template<class T> inline size_t memoryof (T * const e) { return sizeof e; }
template<class T> inline size_t memoryof (T const * const e) { return sizeof e; }
template<class T> inline size_t memoryof (typename list<T>::iterator const & i) { return sizeof i; }
template<class T> inline size_t memoryof (typename list<T>::const_iterator const & i) { return sizeof i; }

template<class T> size_t memoryof (list<T> const & c) CCTK_ATTRIBUTE_PURE;
template<class T> size_t memoryof (set<T> const & c) CCTK_ATTRIBUTE_PURE;
template<class T> size_t memoryof (stack<T> const & c) CCTK_ATTRIBUTE_PURE;
template<class T> size_t memoryof (vector<T> const & c) CCTK_ATTRIBUTE_PURE;



// Container input
template<class T> istream& input (istream& is, list<T>& l);
template<class T> istream& input (istream& is, set<T>& s);
template<class T> istream& input (istream& is, vector<T>& v);

template<class T>
inline istream& operator>> (istream& is, list<T>& l) {
  return input(is,l);
}

template<class T>
inline istream& operator>> (istream& is, set<T>& s) {
  return input(is,s);
}

template<class T>
inline istream& operator>> (istream& is, vector<T>& v) {
  return input(is,v);
}



// Container output
template<class T> ostream& output (ostream& os, const list<T>& l);
template<class S, class T> ostream& output (ostream& os, const map<S,T>& m);
template<class S, class T> ostream& output (ostream& os, const pair<S,T>& p);
template<class T> ostream& output (ostream& os, const set<T>& s);
template<class T> ostream& output (ostream& os, const stack<T>& s);
template<class T> ostream& output (ostream& os, const vector<T>& v);

template<class T>
inline ostream& operator<< (ostream& os, const list<T>& l) {
  return output(os,l);
}

template<class S, class T>
inline ostream& operator<< (ostream& os, const map<S,T>& m) {
  return output(os,m);
}

template<class T>
inline ostream& operator<< (ostream& os, const set<T>& s) {
  return output(os,s);
}

template<class T>
inline ostream& operator<< (ostream& os, const stack<T>& s) {
  return output(os,s);
}

template<class T>
inline ostream& operator<< (ostream& os, const vector<T>& v) {
  return output(os,v);
}



#endif // DEFS_HH
