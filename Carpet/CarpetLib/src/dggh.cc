/***************************************************************************
                          dggh.cc  -  Dimension Generic Grid Hierarchy
                          bounding boxes for each multigrid level of each
                          component of each refinement level
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/Attic/dggh.cc,v 1.2 2002/05/05 22:17:01 schnetter Exp $

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
#include <iostream>

#include "defs.hh"
#include "dh.hh"
#include "th.hh"

#include "dggh.hh"

using namespace std;



  // Constructors
dimgeneric_gh::dimgeneric_gh (const int reffact, const centering refcent,
			      const int mgfact, const centering mgcent)
  : reffact(reffact), refcent(refcent),
    mgfact(mgfact), mgcent(mgcent)
{
  assert (reffact>=1);
  assert (mgfact>=1);
  assert (refcent==vertex_centered || refcent==cell_centered);
  assert (mgcent==vertex_centered || mgcent==cell_centered);
}

// Destructors
dimgeneric_gh::~dimgeneric_gh () { }

// Time hierarchy management
void dimgeneric_gh::add (th* t) {
  ths.push_back(t);
}

void dimgeneric_gh::remove (th* t) {
  ths.remove(t);
}
