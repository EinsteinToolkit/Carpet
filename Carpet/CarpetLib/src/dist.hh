#ifndef DIST_HH
#define DIST_HH

#include <cassert>
#include <cstdio>
#include <cstdlib>

#if 0
#include <complex>
#endif

#include <mpi.h>

#include "cctk.h"

#include "defs.hh"

using namespace std;



namespace dist {
  
  extern MPI_Comm comm;
  
#if 0
  extern MPI_Datatype mpi_complex_float;
  extern MPI_Datatype mpi_complex_double;
  extern MPI_Datatype mpi_complex_long_double;
#else
  extern MPI_Datatype mpi_complex8;
  extern MPI_Datatype mpi_complex16;
  extern MPI_Datatype mpi_complex32;
#endif
  
  void init (int& argc, char**& argv);
  void pseudoinit ();
  void finalize ();
  
  // Debugging output
#define CHECKPOINT dist::checkpoint(__FILE__, __LINE__)
  void checkpoint (const char* file, int line);
  
  
  
  // Information about the communicator
  
  // Rank in the communicator (this processor's number, 0 .. size-1)
  inline int rank ()
  {
    int rank_;
    MPI_Comm_rank (comm, &rank_);
    return rank_;
  }
  
  // Size of the communicator
  inline int size ()
  {
    int size_;
    MPI_Comm_size (comm, &size_);
    return size_;
  }
  
  
  
  // Datatype helpers
  inline MPI_Datatype datatype (const char&)
  { return MPI_CHAR; }
  
  inline MPI_Datatype datatype (const signed char&)
  { return MPI_UNSIGNED_CHAR; }
  
  inline MPI_Datatype datatype (const unsigned char&)
  { return MPI_BYTE; }
  
  inline MPI_Datatype datatype (const short&)
  { return MPI_SHORT; }
  
  inline MPI_Datatype datatype (const unsigned short&)
  { return MPI_UNSIGNED_SHORT; }
  
  inline MPI_Datatype datatype (const int&)
  { return MPI_INT; }
  
  inline MPI_Datatype datatype (const unsigned int&)
  { return MPI_UNSIGNED; }
  
  inline MPI_Datatype datatype (const long&)
  { return MPI_LONG; }
  
  inline MPI_Datatype datatype (const unsigned long&)
  { return MPI_UNSIGNED_LONG; }
  
  inline MPI_Datatype datatype (const long long&)
  { return MPI_LONG_LONG_INT; }
  
  inline MPI_Datatype datatype (const float&)
  { return MPI_FLOAT; }
  
  inline MPI_Datatype datatype (const double&)
  { return MPI_DOUBLE; }
  
  inline MPI_Datatype datatype (const long double&)
  { return MPI_LONG_DOUBLE; }
  
#if 0
  
  inline MPI_Datatype datatype (const complex<float>&)
  { return mpi_complex_float; }
  
  inline MPI_Datatype datatype (const complex<double>&)
  { return mpi_complex_double; }
  
  inline MPI_Datatype datatype (const complex<long double>&)
  { return mpi_complex_long_double; }
  
#else
  
#  ifdef CCTK_REAL4
  inline MPI_Datatype datatype (const CCTK_COMPLEX8&)
  { return mpi_complex8; }
#  endif
  
#  ifdef CCTK_REAL8
  inline MPI_Datatype datatype (const CCTK_COMPLEX16&)
  { return mpi_complex16; }
#  endif
  
#  ifdef CCTK_REAL16
  inline MPI_Datatype datatype (const CCTK_COMPLEX32&)
  { return mpi_complex32; }
#  endif
  
#endif
  
} // namespace dist



#endif // DIST_HH
