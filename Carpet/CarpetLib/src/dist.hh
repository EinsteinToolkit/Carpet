/***************************************************************************
                          dist.hh  -  Helpers for distributed computing
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/dist.hh,v 1.5 2002/05/05 22:17:01 schnetter Exp $

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DIST_HH
#define DIST_HH

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <complex>

#include <mpi.h>

#include "defs.hh"

using namespace std;



namespace dist {
  
  const int tag = 1;
  
  extern MPI_Comm comm;
  
  extern MPI_Datatype mpi_complex_float;
  extern MPI_Datatype mpi_complex_double;
  extern MPI_Datatype mpi_complex_long_double;
  
  void init (int& argc, char**& argv);
  void pseudoinit ();
  void finalize ();
  
  // Debugging output
#define CHECKPOINT dist::checkpoint(__FILE__, __LINE__)
  void checkpoint (const char* file, int line);
  
  
  
  // Datatype helpers
  
  // This generic routine is only declared and not defined.  Only
  // specialised versions have definitions.  Type errors will be
  // caught by the compiler.
  template<class T>
  MPI_Datatype datatype (const T& dummy);
  
  template<class T>
  inline MPI_Datatype datatype (const T& dummy)
  { abort(); return -1; }
  
  template<>
  inline MPI_Datatype datatype (const char& dummy)
  { return MPI_CHAR; }
  
  template<>
  inline MPI_Datatype datatype (const signed char& dummy)
  { return MPI_UNSIGNED_CHAR; }
  
  template<>
  inline MPI_Datatype datatype (const unsigned char& dummy)
  { return MPI_BYTE; }
  
  template<>
  inline MPI_Datatype datatype (const short& dummy)
  { return MPI_SHORT; }
  
  template<>
  inline MPI_Datatype datatype (const unsigned short& dummy)
  { return MPI_UNSIGNED_SHORT; }
  
  template<>
  inline MPI_Datatype datatype (const int& dummy)
  { return MPI_INT; }
  
  template<>
  inline MPI_Datatype datatype (const unsigned int& dummy)
  { return MPI_UNSIGNED; }
  
  template<>
  inline MPI_Datatype datatype (const long& dummy)
  { return MPI_LONG; }
  
  template<>
  inline MPI_Datatype datatype (const unsigned long& dummy)
  { return MPI_UNSIGNED_LONG; }
  
  template<>
  inline MPI_Datatype datatype (const long long& dummy)
  { return MPI_LONG_LONG_INT; }
  
  template<>
  inline MPI_Datatype datatype (const float& dummy)
  { return MPI_FLOAT; }
  
  template<>
  inline MPI_Datatype datatype (const double& dummy)
  { return MPI_DOUBLE; }
  
  template<>
  inline MPI_Datatype datatype (const long double& dummy)
  { return MPI_LONG_DOUBLE; }
  
  template<>
  inline MPI_Datatype datatype (const complex<float>& dummy)
  { return mpi_complex_float; }
  
  template<>
  inline MPI_Datatype datatype (const complex<double>& dummy)
  { return mpi_complex_double; }
  
  template<>
  inline MPI_Datatype datatype (const complex<long double>& dummy)
  { return mpi_complex_long_double; }
  
} // namespace dist



#endif // DIST_HH
