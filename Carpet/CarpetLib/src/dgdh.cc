/***************************************************************************
                          dgdh.cc  -  Dimension Generic Data Hierarchy
                          A grid hierarchy plus ghost zones
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/Attic/dgdh.cc,v 1.4 2002/10/24 11:36:34 schnetter Exp $

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

#include "dgdh.hh"

using namespace std;



// Constructors
dimgeneric_dh::dimgeneric_dh (int prolongation_order_space)
  : prolongation_order_space(prolongation_order_space)
{
  assert (prolongation_order_space>=0);
}

// Destructors
dimgeneric_dh::~dimgeneric_dh ()
{
}

// Helpers
int dimgeneric_dh::prolongation_stencil_size () const {
  assert (prolongation_order_space>=0);
  return prolongation_order_space/2;
}
