/***************************************************************************
                          vect.cc  -  Small inline vectors
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/vect.cc,v 1.4 2001/03/27 22:26:31 eschnett Exp $

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
template<class T,int D>
void vect<T,D>::output (ostream& os) const {
  os << "[";
  for (int d=0; d<D; ++d) {
    if (d>0) os << ",";
    os << (*this)[d];
  }
  os << "]";
}



#if defined(TMPL_EXPLICIT)
template class vect<int,1>;
template class vect<int,2>;
template class vect<int,3>;
#endif
