/***************************************************************************
                          gf.cc  -  Grid Function
                          data for every element of a data hierarchy
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/gf.cc,v 1.3 2001/03/22 18:42:05 eschnett Exp $

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

#include "defs.hh"

#if !defined(TMPL_IMPLICIT) || !defined(GF_HH)
#  include "gf.hh"
#endif

using namespace std;



// Constructors
template<class T,int D>
gf<T,D>::gf (const string name, th<D>& t, dh<D>& d,
	     const int tmin, const int tmax)
  : generic_gf<D>(name, t, d, tmin, tmax)
{
  recompose();
}

// Destructors
template<class T,int D>
gf<T,D>::~gf () { }



// Access to the data
template<class T,int D>
const data<T,D>* gf<T,D>::operator() (int tl, int rl, int c, int ml) const {
  assert (tl>=tmin && tl<=tmax);
  assert (rl>=0 && rl<h.reflevels());
  assert (c>=0 && c<h.components(rl));
  assert (ml>=0 && ml<h.mglevels(rl,c));
  return (const data<T,D>*)storage[tl-tmin][rl][c][ml];
}

template<class T,int D>
data<T,D>* gf<T,D>::operator() (int tl, int rl, int c, int ml) {
  assert (tl>=tmin && tl<=tmax);
  assert (rl>=0 && rl<h.reflevels());
  assert (c>=0 && c<h.components(rl));
  assert (ml>=0 && ml<h.mglevels(rl,c));
  return (data<T,D>*)storage[tl-tmin][rl][c][ml];
}



// Output
template<class T,int D>
ostream& gf<T,D>::out (ostream& os) const {
  os << "gf<" STR(T) "," << D << ">:\"" << name << "\","
     << "dt=[" << tmin << ":" << tmax<< "]";
  return os;
}



#if defined(TMPL_EXPLICIT)

#define INSTANTIATE(T)				\
template class gf<T,3>;

#include "instantiate"

#undef INSTANTIATE

#endif
