/***************************************************************************
                          vect.cc  -  Small inline vectors
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/vect.cc,v 1.3 2001/03/22 18:42:06 eschnett Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <assert.h>

#include <iostream>

#include "defs.hh"

#if !defined(TMPL_IMPLICIT) || !defined(VECT_HH)
#  include "vect.hh"
#endif

using namespace std;



// Output
#ifndef SGI
// This doesn't work on SGIs.  Is this legal C++?
template<class T,int D>
ostream& operator<< (ostream& os, const vect<T,D>& a) {
  os << "[";
  for (int d=0; d<D; ++d) {
    if (d>0) os << ",";
    os << a[d];
  }
  os << "]";
  return os;
}
#else
ostream& operator<< (ostream& os, const vect<int,1>& a) {
  os << "[";
  for (int d=0; d<1; ++d) {
    if (d>0) os << ",";
    os << a[d];
  }
  os << "]";
  return os;
}
ostream& operator<< (ostream& os, const vect<int,2>& a) {
  os << "[";
  for (int d=0; d<2; ++d) {
    if (d>0) os << ",";
    os << a[d];
  }
  os << "]";
  return os;
}
ostream& operator<< (ostream& os, const vect<int,3>& a) {
  os << "[";
  for (int d=0; d<3; ++d) {
    if (d>0) os << ",";
    os << a[d];
  }
  os << "]";
  return os;
}
#endif



#if defined(TMPL_EXPLICIT)
template class vect<int,1>;
template ostream& operator<< (ostream& os, const vect<int,1>& a);

template class vect<int,2>;
template ostream& operator<< (ostream& os, const vect<int,2>& a);

template class vect<int,3>;
template ostream& operator<< (ostream& os, const vect<int,3>& a);
#endif
