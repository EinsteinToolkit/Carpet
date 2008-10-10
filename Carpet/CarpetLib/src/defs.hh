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
#include <set>
#include <stack>
#include <vector>

#include "cctk.h"



using namespace std;


  
// A compile time pseudo assert statement
#define static_assert(_x, _msg) do { typedef int ai[(_x) ? 1 : -1]; } while(0)



// Check a return value
#define check(_expr) do { bool const _val = (_expr); assert(_val); } while(0)



// Define the restrict qualifier
#ifdef CCTK_CXX_RESTRICT
#  define restrict CCTK_CXX_RESTRICT
#endif



// Use this macro AT instead of vector's operator[] or at().
// Depending on the macro NDEBUG, this macro AT either checks for
// valid indices or not.
#ifndef CARPET_OPTIMISE
#  define AT(index) at(index)
#else
#  define AT(index) operator[](index)
#endif



// Begin a new line without flushing the output buffer
char const * const eol = "\n";



// Number of dimensions
const int dim = 3;
  


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
inline T square (const T& x) { return x*x; }

// Another useful helper
template<class T>
T ipow (T x, int y);



// Input streams
struct input_error { };
void skipws (istream& is);
void expect (istream& is, char c);
void consume (istream& is, char c);
void consume (istream& is, char const * c);



// Names for types

#ifdef HAVE_CCTK_INT1
inline const char * typestring (const CCTK_INT1& dummy)
{ return "CCTK_INT1"; }
#endif

#ifdef HAVE_CCTK_INT2
inline const char * typestring (const CCTK_INT2& dummy)
{ return "CCTK_INT2"; }
#endif

#ifdef HAVE_CCTK_INT4
inline const char * typestring (const CCTK_INT4& dummy)
{ return "CCTK_INT4"; }
#endif

#ifdef HAVE_CCTK_INT8
inline const char * typestring (const CCTK_INT8& dummy)
{ return "CCTK_INT8"; }
#endif

#ifdef HAVE_CCTK_REAL4
inline const char * typestring (const CCTK_REAL4& dummy)
{ return "CCTK_REAL4"; }
#endif

#ifdef HAVE_CCTK_REAL8
inline const char * typestring (const CCTK_REAL8& dummy)
{ return "CCTK_REAL8"; }
#endif

#ifdef HAVE_CCTK_REAL16
inline const char * typestring (const CCTK_REAL16& dummy)
{ return "CCTK_REAL16"; }
#endif

#ifdef HAVE_CCTK_REAL4
inline const char * typestring (const CCTK_COMPLEX8& dummy)
{ return "CCTK_COMPLEX8"; }
#endif

#ifdef HAVE_CCTK_REAL8
inline const char * typestring (const CCTK_COMPLEX16& dummy)
{ return "CCTK_COMPLEX16"; }
#endif

#ifdef HAVE_CCTK_REAL16
inline const char * typestring (const CCTK_COMPLEX32& dummy)
{ return "CCTK_COMPLEX32"; }
#endif



namespace CarpetLib {
  namespace good {
    
    // Explicitly overload abs for all types in the same namespace, to
    // circumvent confusion among some compilers
    
    // CCTK_BYTE is unsigned
    inline CCTK_BYTE abs (CCTK_BYTE const & x) { return x; }
    
#if 0
    // This does not work on AIX, which does not have long long abs
    // (long long)
#  ifdef HAVE_CCTK_INT1
    inline CCTK_INT1 abs (CCTK_INT1 const & x) { return std::abs (x); }
#  endif
#  ifdef HAVE_CCTK_INT2
    inline CCTK_INT2 abs (CCTK_INT2 const & x) { return std::abs (x); }
#  endif
#  ifdef HAVE_CCTK_INT4
    inline CCTK_INT4 abs (CCTK_INT4 const & x) { return std::abs (x); }
#  endif
#  ifdef HAVE_CCTK_INT8
    inline CCTK_INT8 abs (CCTK_INT8 const & x) { return std::abs (x); }
#  endif
#endif
    
#if 0
    // This does not work on Linux with Intel compilers, which do not
    // always have long long llabs (long long)
    inline signed char abs (signed char const & x) { return ::abs (x); }
    inline unsigned char abs (unsigned char const & x) { return ::abs (x); }
    inline short abs (short const & x) { return ::abs (x); }
    inline int abs (int const & x) { return ::abs (x); }
    inline long abs (long const & x) { return ::labs (x); }
#  ifdef SIZEOF_LONG_LONG
    inline long long abs (long long const & x) { return ::llabs (x); }
#  endif
#endif
    
#if 1
#  ifdef HAVE_CCTK_INT1
    inline CCTK_INT1 abs (CCTK_INT1 const & x) { return x < 0 ? - x : x; }
#  endif
#  ifdef HAVE_CCTK_INT2
    inline CCTK_INT2 abs (CCTK_INT2 const & x) { return x < 0 ? - x : x; }
#  endif
#  ifdef HAVE_CCTK_INT4
    inline CCTK_INT4 abs (CCTK_INT4 const & x) { return x < 0 ? - x : x; }
#  endif
#  ifdef HAVE_CCTK_INT8
    inline CCTK_INT8 abs (CCTK_INT8 const & x) { return x < 0 ? - x : x; }
#  endif
#endif
    
#ifdef HAVE_CCTK_REAL4
    inline CCTK_REAL4 abs (CCTK_REAL4 const & x) { return std::abs (x); }
#endif
#ifdef HAVE_CCTK_REAL8
    inline CCTK_REAL8 abs (CCTK_REAL8 const & x) { return std::abs (x); }
#endif
#ifdef HAVE_CCTK_REAL16
    inline CCTK_REAL16 abs (CCTK_REAL16 const & x) { return std::abs (x); }
#endif
    
#ifdef HAVE_CCTK_COMPLEX8
    inline CCTK_REAL4 abs (CCTK_COMPLEX8 const & x)
    { return CCTK_Cmplx8Abs (x); }
#endif
#ifdef HAVE_CCTK_COMPLEX16
    inline CCTK_REAL8 abs (CCTK_COMPLEX16 const & x)
    { return CCTK_Cmplx16Abs (x); }
#endif
#ifdef HAVE_CCTK_COMPLEX32
    inline CCTK_REAL16 abs (CCTK_COMPLEX32 const & x)
    { return CCTK_Cmplx32Abs (x); }
#endif
    
  } // namespace good
} // namespace CarpetLib



// Container memory usage
inline size_t memoryof (char e) { return sizeof e; }
inline size_t memoryof (short e) { return sizeof e; }
inline size_t memoryof (int e) { return sizeof e; }
inline size_t memoryof (long e) { return sizeof e; }
inline size_t memoryof (long long e) { return sizeof e; }
inline size_t memoryof (unsigned char e) { return sizeof e; }
inline size_t memoryof (unsigned short e) { return sizeof e; }
inline size_t memoryof (unsigned int e) { return sizeof e; }
inline size_t memoryof (unsigned long e) { return sizeof e; }
inline size_t memoryof (unsigned long long e) { return sizeof e; }
inline size_t memoryof (float e) { return sizeof e; }
inline size_t memoryof (double e) { return sizeof e; }
inline size_t memoryof (long double e) { return sizeof e; }
inline size_t memoryof (void * e) { return sizeof e; }
template<class T> inline size_t memoryof (T * e) { return sizeof e; }
template<class T> inline size_t memoryof (T const * e) { return sizeof e; }
template<class T> size_t memoryof (list<T> const & c);
template<class T> size_t memoryof (set<T> const & c);
template<class T> size_t memoryof (stack<T> const & c);
template<class T> size_t memoryof (vector<T> const & c);



// Container input
template<class T> istream& input (istream& is, vector<T>& v);

template<class T>
inline istream& operator>> (istream& is, vector<T>& v) {
  return input(is,v);
}



// Container output
template<class T> ostream& output (ostream& os, const list<T>& l);
template<class T> ostream& output (ostream& os, const set<T>& s);
template<class T> ostream& output (ostream& os, const stack<T>& s);
template<class T> ostream& output (ostream& os, const vector<T>& v);

template<class T>
inline ostream& operator<< (ostream& os, const list<T>& l) {
  return output(os,l);
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
