/***************************************************************************
                          dgdata.cc  -  description
                             -------------------
    begin                : Wed Jul 19 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/Attic/dgdata.cc,v 1.2 2002/05/05 22:17:00 schnetter Exp $

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

#include "dgdata.hh"

using namespace std;



// Constructors
dimgeneric_data::dimgeneric_data ()
  : _has_storage(false)
{ }

// Destructors
dimgeneric_data::~dimgeneric_data () { }
