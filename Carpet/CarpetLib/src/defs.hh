/***************************************************************************
                          defs.hh  -  Commonly used definitions
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/defs.hh,v 1.2 2001/03/07 13:00:57 eschnett Exp $

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

#include <algorithm>
#include <list>
#include <set>
#include <vector>

// Stringification
#define STR(s) #s

// Fortranification
#define FORTRAN_NAME(x) x##_

// Fortran style function arguments
#define restrict __restrict__

// A general type
enum centering { vertex_centered, cell_centered };

// Useful helper
template<class T>
inline T square (const T& x) { return x*x; }

// Container output
template<class T> ostream& operator<< (ostream& os, const list<T>& l);
template<class T> ostream& operator<< (ostream& os, const set<T>& s);
template<class T> ostream& operator<< (ostream& os, const vector<T>& v);



#if defined(TMPL_IMPLICIT)
#  include "defs.cc"
#endif

#endif // DEFS_HH
