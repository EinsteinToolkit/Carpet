/***************************************************************************
                          dggf.cc  -  Dimension Generic Grid Function
                          grid function without type information
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/Attic/dggf.cc,v 1.1 2001/06/12 14:56:58 schnetter Exp $

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
#include <stdlib.h>

// #include <iostream>
// #include <string>

// #include "defs.hh"
// #include "dh.hh"
// #include "th.hh"

#if !defined(TMPL_IMPLICIT) || !defined(GGF_HH)
#  include "dggf.hh"
#endif

using namespace std;



// Constructors
dimgeneric_gf::dimgeneric_gf (const string name, th& t,
			      const int tmin, const int tmax)
  : name(name), t(t), tmin(tmin), tmax(tmax)
{
  assert (tmin<=tmax+1);
}

// Destructors
dimgeneric_gf::~dimgeneric_gf ()
{
}

// Comparison
bool dimgeneric_gf::operator== (const dimgeneric_gf& f) const {
  return this == &f;
}



#if defined(TMPL_EXPLICIT)
#endif
