#ifndef DEFS_HH
#define DEFS_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>
#include <cassert>
#include <iostream>
#include <list>
#include <set>
#include <stack>
#include <vector>

#include "cctk.h"



using namespace std;


  
// A compile time pseudo assert statement
#define static_assert(_x) do { typedef int ai[(_x) ? 1 : -1]; } while(0)



// Define the restrict qualifier
#ifdef CCTK_CXX_RESTRICT
#  define restrict CCTK_CXX_RESTRICT
#endif



// Use this macro AT instead of vector's operator[] or at().
// Depending on the macro NDEBUG, this macro AT either checks for
// valid indices or not.
#ifndef NDEBUG
#  define AT(index) at(index)
#else
#  define AT(index) operator[](index)
#endif



// Begin a new line without flushing the output buffer
char const * const eol = "\n";



// Number of dimensions
const int dim = 3;
  


// Some shortcuts for type names
template<typename T, int D> class bbox;
template<typename T, int D> class bboxset;
template<typename T, int D> class vect;

typedef vect<bool,dim>   bvect;
typedef vect<int,dim>    ivect;
typedef bbox<int,dim>    ibbox;
typedef bboxset<int,dim> ibset;

typedef vect<vect<bool,2>,dim> bbvect;
typedef vect<vect<int,2>,dim>  iivect;

typedef vect<vect<bool,dim>,2> b2vect;
typedef vect<vect<int,dim>,2>  i2vect;



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
