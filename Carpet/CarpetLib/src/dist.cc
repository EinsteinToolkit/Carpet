/***************************************************************************
                          dist.cc  -  Helpers for distributed computing
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/dist.cc,v 1.1 2001/03/01 13:40:10 eschnett Exp $

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

#include <mpi.h>

#include "defs.hh"

#if !defined(TMPL_IMPLICIT) || !defined(DIST_HH)
#  include "dist.hh"
#endif



MPI_Comm dist::comm;



void dist::init (int& argc, char**& argv) {
  MPI_Init (&argc, &argv);
  comm = MPI_COMM_WORLD;
}

void dist::pseudoinit () {
  comm = MPI_COMM_WORLD;
}

void dist::finalize () {
  MPI_Finalize ();
}



#if defined(TMPL_EXPLICIT)
#endif
