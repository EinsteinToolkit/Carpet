/***************************************************************************
                          dgdh.cc  -  Dimension Generic Data Hierarchy
                          A grid hierarchy plus ghost zones
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/Attic/dgdh.cc,v 1.1 2001/06/12 14:56:57 schnetter Exp $

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
#include "dist.hh"
#include "ggf.hh"

#if !defined(TMPL_IMPLICIT) || !defined(DH_HH)
#  include "dgdh.hh"
#endif

#undef DEBUG_OUTPUT

using namespace std;



// Constructors
dimgeneric_dh::dimgeneric_dh (int prolongation_order)
  : prolongation_order(prolongation_order)
{
  assert (prolongation_order>=0);
}

// Destructors
dimgeneric_dh::~dimgeneric_dh ()
{
}

// Helpers
int dimgeneric_dh::prolongation_stencil_size () const {
  assert (prolongation_order>=0);
  return prolongation_order/2;
}



#if defined(TMPL_EXPLICIT)
#endif
