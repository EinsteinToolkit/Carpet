/***************************************************************************
                          vect.cc  -  Small inline vectors
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/vect.cc,v 1.8 2002/05/05 22:17:03 schnetter Exp $

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

#include "vect.hh"

using namespace std;



// Input
template<class T,int D>
void vect<T,D>::input (istream& is) {
  skipws (is);
  assert (is.peek() == '[');
  is.get();
  for (int d=0; d<D; ++d) {
    is >> (*this)[d];
    if (d<D-1) {
      skipws (is);
      assert (is.peek() == ',');
      is.get();
    }
  }
  skipws (is);
  assert (is.peek() == ']');
  is.get();
}



// Output
template<class T,int D>
void vect<T,D>::output (ostream& os) const {
  os << "[";
  for (int d=0; d<D; ++d) {
    os << (*this)[d];
    if (d<D-1) os << ",";
  }
  os << "]";
}



// Note: We need all dimensions all the time.
template class vect<int,1>;
template class vect<int,2>;
template class vect<int,3>;

template void vect<double,3>::input (istream& is);
template void vect<vect<bool,2>,3>::input (istream& is);
template void vect<double,3>::output (ostream& os) const;
template void vect<vect<bool,2>,3>::output (ostream& os) const;
