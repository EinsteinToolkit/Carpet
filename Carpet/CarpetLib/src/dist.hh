// $Header: /home/eschnett/C/carpet/Carpet/Carpet/CarpetLib/src/dist.hh,v 1.8 2003/01/03 15:49:36 schnetter Exp $

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
  inline MPI_Datatype datatype (const char& dummy)
  { return MPI_CHAR; }
  
  inline MPI_Datatype datatype (const signed char& dummy)
  { return MPI_UNSIGNED_CHAR; }
  
  inline MPI_Datatype datatype (const unsigned char& dummy)
  { return MPI_BYTE; }
  
  inline MPI_Datatype datatype (const short& dummy)
  { return MPI_SHORT; }
  
  inline MPI_Datatype datatype (const unsigned short& dummy)
  { return MPI_UNSIGNED_SHORT; }
  
  inline MPI_Datatype datatype (const int& dummy)
  { return MPI_INT; }
  
  inline MPI_Datatype datatype (const unsigned int& dummy)
  { return MPI_UNSIGNED; }
  
  inline MPI_Datatype datatype (const long& dummy)
  { return MPI_LONG; }
  
  inline MPI_Datatype datatype (const unsigned long& dummy)
  { return MPI_UNSIGNED_LONG; }
  
  inline MPI_Datatype datatype (const long long& dummy)
  { return MPI_LONG_LONG_INT; }
  
  inline MPI_Datatype datatype (const float& dummy)
  { return MPI_FLOAT; }
  
  inline MPI_Datatype datatype (const double& dummy)
  { return MPI_DOUBLE; }
  
  inline MPI_Datatype datatype (const long double& dummy)
  { return MPI_LONG_DOUBLE; }
  
  inline MPI_Datatype datatype (const complex<float>& dummy)
  { return mpi_complex_float; }
  
  inline MPI_Datatype datatype (const complex<double>& dummy)
  { return mpi_complex_double; }
  
  inline MPI_Datatype datatype (const complex<long double>& dummy)
  { return mpi_complex_long_double; }
  
} // namespace dist



#endif // DIST_HH
