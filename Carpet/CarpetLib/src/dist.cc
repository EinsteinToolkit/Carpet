/***************************************************************************
                          dist.cc  -  Helpers for distributed computing
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/dist.cc,v 1.3 2001/03/07 13:00:57 eschnett Exp $

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

#include "cctk.h"
#include "cctk_Parameters.h"

#include "defs.hh"

#if !defined(TMPL_IMPLICIT) || !defined(DIST_HH)
#  include "dist.hh"
#endif



namespace dist {
  
  MPI_Comm comm;
  
  MPI_Datatype mpi_complex_float;
  MPI_Datatype mpi_complex_double;
  MPI_Datatype mpi_complex_long_double;
  
  void init (int& argc, char**& argv) {
    MPI_Init (&argc, &argv);
    pseudoinit();
  }
  
  void pseudoinit () {
    comm = MPI_COMM_WORLD;
    
    MPI_Type_contiguous (2, MPI_FLOAT, &mpi_complex_float);
    MPI_Type_commit (&mpi_complex_float);
    MPI_Type_contiguous (2, MPI_DOUBLE, &mpi_complex_double);
    MPI_Type_commit (&mpi_complex_double);
    MPI_Type_contiguous (2, MPI_LONG_DOUBLE, &mpi_complex_long_double);
    MPI_Type_commit (&mpi_complex_long_double);
  }
  
  void finalize () {
    MPI_Finalize ();
  }
  
  void checkpoint (const char* file, int line) {
    DECLARE_CCTK_PARAMETERS;
    if (verbose) {
      int rank;
      MPI_Comm_rank (comm, &rank);
      printf ("CHECKPOINT: processor %d, file %s, line %d\n",
	      rank, file, line);
    }
    if (barriers) {
      MPI_Barrier (comm);
    }
  }
  
} // namespace dist



#if defined(TMPL_EXPLICIT)
#endif
