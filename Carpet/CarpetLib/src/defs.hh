/***************************************************************************
                          defs.hh  -  Commonly used definitions
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/defs.hh,v 1.6 2002/01/09 13:56:26 schnetter Exp $

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



#if defined(TMPL_IMPLICIT)
#  include "defs.cc"
#endif

#endif // DEFS_HH
