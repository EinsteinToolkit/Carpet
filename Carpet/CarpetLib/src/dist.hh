/***************************************************************************
                          dist.hh  -  Helpers for distributed computing
                             -------------------
    begin                : Sun Jun 11 2000
    copyright            : (C) 2000 by Erik Schnetter
    email                : schnetter@astro.psu.edu

    $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/dist.hh,v 1.2 2001/03/05 21:48:38 eschnett Exp $

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

#include <cassert>
#include <cstdio>
#include <cstdlib>

#include <mpi.h>

#include "defs.hh"



// A checkpoint for debugging purposes
#define DIST_VERBOSE						\
do {								\
  int rank;							\
  MPI_Comm_rank (dist::comm, &rank);				\
  printf ("CHECKPOINT: processor %d, file %s, line %d\n",	\
	  rank, __FILE__, __LINE__);				\
} while(0)

// A barrier for debugging purposes
#define DIST_BARRIER				\
do {						\
  MPI_Barrier(dist::comm);			\
} while(0)

// Both of the above
#define DIST_VERBOSE_BARRIER			\
do {						\
  DIST_VERBOSE;					\
  DIST_BARRIER;					\
} while(0)

// Do nothing
#define DIST_NODEBUG do {} while(0)



struct dist {
  
  static const int tag = 1;
  
  static MPI_Comm comm;
  
  static void init (int& argc, char**& argv);
  static void pseudoinit ();
  static void finalize ();
  
  template<class T>
  static MPI_Datatype datatype (const T& dummy);
  
};



template<class T>
inline MPI_Datatype dist::datatype (const T& dummy)
{ abort(); return -1; }

template<>
inline MPI_Datatype dist::datatype (const char& dummy)
{ return MPI_CHAR; }

template<>
inline MPI_Datatype dist::datatype (const signed char& dummy)
{ return MPI_UNSIGNED_CHAR; }
  
template<>
inline MPI_Datatype dist::datatype (const unsigned char& dummy)
{ return MPI_BYTE; }

template<>
inline MPI_Datatype dist::datatype (const short& dummy)
{ return MPI_SHORT; }

template<>
inline MPI_Datatype dist::datatype (const unsigned short& dummy)
{ return MPI_UNSIGNED_SHORT; }

template<>
inline MPI_Datatype dist::datatype (const int& dummy)
{ return MPI_INT; }

template<>
inline MPI_Datatype dist::datatype (const unsigned int& dummy)
{ return MPI_UNSIGNED; }

template<>
inline MPI_Datatype dist::datatype (const long& dummy)
{ return MPI_LONG; }

template<>
inline MPI_Datatype dist::datatype (const unsigned long& dummy)
{ return MPI_UNSIGNED_LONG; }
  
template<>
inline MPI_Datatype dist::datatype (const long long& dummy)
{ return MPI_LONG_LONG_INT; }

template<>
inline MPI_Datatype dist::datatype (const float& dummy)
{ return MPI_FLOAT; }

template<>
inline MPI_Datatype dist::datatype (const double& dummy)
{ return MPI_DOUBLE; }

template<>
inline MPI_Datatype dist::datatype (const long double& dummy)
{ return MPI_LONG_DOUBLE; }



#if defined(TMPL_IMPLICIT)
#  include "dist.cc"
#endif

#endif // DIST_HH
