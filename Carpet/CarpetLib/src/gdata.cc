/***************************************************************************
                          gdata.cc  -  description
                             -------------------
    begin                : Wed Jul 19 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/gdata.cc,v 1.1 2001/03/01 13:40:10 eschnett Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <cassert>
#include <fstream>

#include "bbox.hh"
#include "defs.hh"
#include "dist.hh"
#include "vect.hh"

#if !defined(TMPL_IMPLICIT) || !defined(GDATA_HH)
#  include "gdata.hh"
#endif



// Constructors
template<int D>
generic_data<D>::generic_data ()
  : _has_storage(false)
{ }

// Destructors
template<int D>
generic_data<D>::~generic_data () { }



// Output
template<int D>
ostream& operator<< (ostream& os, const generic_data<D>& f) {
  return f.out(os);
}



#if defined(TMPL_EXPLICIT)
template class generic_data<1>;
template ostream& operator<< (ostream& os, const generic_data<1>& d);

template class generic_data<2>;
template ostream& operator<< (ostream& os, const generic_data<2>& d);

template class generic_data<3>;
template ostream& operator<< (ostream& os, const generic_data<3>& d);
#endif
