/***************************************************************************
                          defs.hh  -  Commonly used definitions
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/defs.hh,v 1.9 2002/12/31 13:29:07 schnetter Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DEFS_HH
#define DEFS_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>

#include <algorithm>
#include <complex>
#include <iostream>
#include <list>
#include <set>
#include <vector>

using namespace std;

// Stringification
#define STR(s) #s

// Fortranification
#define FORTRAN_NAME(x) x##_

// A general type
enum centering { vertex_centered, cell_centered };

// Useful helper
template<class T>
inline T square (const T& x) { return x*x; }

// Another useful helper
template<class T>
inline T ipow (const T& x, const int y) {
  if (y<0) {
    return T(1)/ipow(x,-y);
  } else if (y==0) {
    return T(1);
  } else if (y%2) {
    return x * ipow(x*x,y/2);
  } else {
    return ipow(x*x,y/2);
  }
}



// Skip whitespace
void skipws (istream& is);



// Names for types
inline const char * typestring (const char& dummy)
{ return "char"; }

inline const char * typestring (const signed char& dummy)
{ return "signed char"; }

inline const char * typestring (const unsigned char& dummy)
{ return "unsigned char"; }

inline const char * typestring (const short& dummy)
{ return "short"; }

inline const char * typestring (const unsigned short& dummy)
{ return "unsigned short"; }

inline const char * typestring (const int& dummy)
{ return "int"; }

inline const char * typestring (const unsigned int& dummy)
{ return "unsigned int"; }

inline const char * typestring (const long& dummy)
{ return "long"; }

inline const char * typestring (const unsigned long& dummy)
{ return "unsigned long"; }

inline const char * typestring (const long long& dummy)
{ return "long long"; }

inline const char * typestring (const unsigned long long& dummy)
{ return "unsigned long long"; }

inline const char * typestring (const float& dummy)
{ return "float"; }

inline const char * typestring (const double& dummy)
{ return "double"; }

inline const char * typestring (const long double& dummy)
{ return "long double"; }

inline const char * typestring (const complex<float>& dummy)
{ return "complex<float>"; }

inline const char * typestring (const complex<double>& dummy)
{ return "complex<double>"; }

inline const char * typestring (const complex<long double>& dummy)
{ return "complex<long double>"; }



// Container input
template<class T> istream& input (istream& is, vector<T>& v);

template<class T>
inline istream& operator>> (istream& is, vector<T>& v) {
  return input(is,v);
}



// Container output
template<class T> ostream& output (ostream& os, const list<T>& l);
template<class T> ostream& output (ostream& os, const set<T>& s);
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
inline ostream& operator<< (ostream& os, const vector<T>& v) {
  return output(os,v);
}



#endif // DEFS_HH
